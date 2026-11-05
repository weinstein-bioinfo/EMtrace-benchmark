"""TSV / JSON / Excel writers for the detection outputs (§2.3)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

from ngs_benchmark.schema import TAXGROUP_SHORT, TAXRANKS


# Columns in the SUMMARY tab that should be formatted as percentages.
_PERCENT_COLS = ("Pr", "Rc", "F1", "mean_accuracy")
# Red fill matches the legacy R fgFill="#ff5050" used to mark FP/FN rows.
_FALSE_FILL = PatternFill(fgColor="FF5050", fill_type="solid")
_HEADER_FONT = Font(bold=True)


def write_tsv(df: pd.DataFrame, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def _json_default(obj):
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None if math.isnan(obj) else ("Infinity" if obj > 0 else "-Infinity")
    if hasattr(obj, "tolist"):
        return obj.tolist()
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Not JSON serialisable: {type(obj).__name__}")


def write_json(obj: dict, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=_json_default, sort_keys=True)


def write_xlsx_summary(
    metrics_df: pd.DataFrame,
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    path: str | Path,
    taxgroups: tuple[str, ...] = ("Viruses", "Bacteria"),
    species_rank: str = "species",
    top_rank: str = "division",
) -> None:
    """Write the legacy-R `summary.xlsx` for a single (pred, truth) run.

    Tab layout (matches `bmark_stats()` in NGSlib.R):

        SUMMARY                — `metrics_df` with percent format on Pr/Rc/F1/mean_accuracy
        {short}_true           — truth rows in `taxgroup`, FN rows red-filled
        {short}_pred           — pred rows in `taxgroup`, FP rows red-filled
        {short}_TP             — pred rows whose species ∈ truth species set
        {short}_FN             — truth rows whose species ∉ pred species set
        {short}_FP             — pred rows whose species ∉ truth species set

    where `short` is the lower-case 2-letter taxgroup short label
    (`vi`, `ba`, `euk`, ...) from `schema.TAXGROUP_SHORT`.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    wb = Workbook()
    wb.remove(wb.active)  # drop the default sheet

    # SUMMARY ----------------------------------------------------------------
    ws = wb.create_sheet("SUMMARY")
    _write_dataframe(ws, metrics_df)
    _format_percent_cols(ws, metrics_df, _PERCENT_COLS)
    _autosize(ws, metrics_df)

    # Per-taxgroup detail tabs ----------------------------------------------
    detail_cols = [c for c in ("readn", *TAXRANKS) if c in pred.columns or c in truth.columns]

    for group in taxgroups:
        short = TAXGROUP_SHORT.get(group, group[:3].lower())

        pred_grp = (
            pred[pred[top_rank] == group][detail_cols].sort_values(
                "readn", ascending=False, kind="mergesort"
            )
            if top_rank in pred.columns
            else pred.iloc[0:0][detail_cols]
        )
        true_grp = (
            truth[truth[top_rank] == group][detail_cols].sort_values(
                "readn", ascending=False, kind="mergesort"
            )
            if top_rank in truth.columns
            else truth.iloc[0:0][detail_cols]
        )

        pred_species = set(pred_grp[species_rank].dropna().unique())
        true_species = set(true_grp[species_rank].dropna().unique())
        tp_species = pred_species & true_species
        fp_species = pred_species - true_species
        fn_species = true_species - pred_species

        pred_grp = pred_grp.reset_index(drop=True)
        true_grp = true_grp.reset_index(drop=True)

        # {short}_true ─ truth rows; rows whose species is a FN are red-filled.
        ws_true = wb.create_sheet(f"{short}_true")
        _write_dataframe(ws_true, true_grp)
        _highlight_rows_where(ws_true, true_grp, species_rank, fn_species)
        _format_int_col(ws_true, true_grp, "readn")
        _autosize(ws_true, true_grp)

        # {short}_pred ─ pred rows; FP rows red-filled.
        ws_pred = wb.create_sheet(f"{short}_pred")
        _write_dataframe(ws_pred, pred_grp)
        _highlight_rows_where(ws_pred, pred_grp, species_rank, fp_species)
        _format_int_col(ws_pred, pred_grp, "readn")
        _autosize(ws_pred, pred_grp)

        # Status-filtered tables — TP, FN, FP.
        for tab, source, species_set in (
            (f"{short}_TP", pred_grp, tp_species),
            (f"{short}_FN", true_grp, fn_species),
            (f"{short}_FP", pred_grp, fp_species),
        ):
            sub = source[source[species_rank].isin(species_set)].reset_index(drop=True)
            ws_sub = wb.create_sheet(tab)
            _write_dataframe(ws_sub, sub)
            _format_int_col(ws_sub, sub, "readn")
            _autosize(ws_sub, sub)

    wb.save(path)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _write_dataframe(ws, df: pd.DataFrame) -> None:
    """Write header + rows. Empty DataFrame writes only the header."""
    ws.append(list(df.columns))
    for cell in ws[1]:
        cell.font = _HEADER_FONT
    for _, row in df.iterrows():
        ws.append([_excel_cell(v) for v in row])


def _excel_cell(value):
    """Convert pandas values to types openpyxl handles without complaint."""
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, float) and math.isinf(value):
        return "Infinity" if value > 0 else "-Infinity"
    if hasattr(value, "item"):
        try:
            return value.item()
        except (ValueError, TypeError):
            return str(value)
    return value


def _format_percent_cols(ws, df: pd.DataFrame, cols: tuple[str, ...]) -> None:
    for col in cols:
        if col not in df.columns:
            continue
        col_idx = list(df.columns).index(col) + 1  # 1-based for openpyxl
        for row_idx in range(2, len(df) + 2):       # 1 is the header
            ws.cell(row=row_idx, column=col_idx).number_format = "0.0%"


def _format_int_col(ws, df: pd.DataFrame, col: str) -> None:
    if col not in df.columns:
        return
    col_idx = list(df.columns).index(col) + 1
    for row_idx in range(2, len(df) + 2):
        ws.cell(row=row_idx, column=col_idx).number_format = "0"


def _highlight_rows_where(ws, df: pd.DataFrame, key_col: str, key_values: set) -> None:
    """Apply red fill to every cell in rows whose `key_col` value is in `key_values`."""
    if key_col not in df.columns or not key_values:
        return
    key_idx = list(df.columns).index(key_col)
    n_cols = len(df.columns)
    for i, value in enumerate(df[key_col].tolist()):
        if value in key_values:
            excel_row = i + 2  # +1 for header, +1 for 1-based indexing
            for j in range(1, n_cols + 1):
                ws.cell(row=excel_row, column=j).fill = _FALSE_FILL
    _ = key_idx  # kept for clarity / future use


def _autosize(ws, df: pd.DataFrame, max_width: int = 40) -> None:
    """Best-effort column auto-sizing. openpyxl has no native auto-fit."""
    for j, col in enumerate(df.columns, start=1):
        header_len = len(str(col))
        sample = df[col].head(50).astype(str)
        max_data = int(sample.str.len().max()) if len(sample) > 0 else 0
        width = min(max(header_len, max_data) + 2, max_width)
        ws.column_dimensions[get_column_letter(j)].width = width
