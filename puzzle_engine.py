# puzzle_engine.py
# -*- coding: utf-8 -*-
"""
מנוע החידות: טעינה מסדרת JSON, סינון לפי ordered/finite,
וברירה רנדומלית של פונקציית חידה מתוך registry.
"""
import json
import random
from typing import Dict, List
from puzzles import puzzle_registry

# נתיב ברירת מחדל לקובץ הסדרות
SERIES_PATH = "series.json"


def load_series(path: str = SERIES_PATH) -> List[Dict]:
    """טוען את הקובץ series.json ומחזיר List[Dict]"""
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get_ordered_series(all_series: List[Dict]) -> List[Dict]:
    """מסנן רק את הסדרות שהן ordered ו-ftine"""
    return [s for s in all_series if s.get("ordered") and s.get("finite")]


def generate_random_puzzle(rng: random.Random, all_series: List[Dict]) -> Dict:
    """
    בוחר פונקציית חידה רנדומלית ו-series_entry מתאים,
    מפעיל אותה עם הפרמטרים הנחוצים.
    """
    # סינון סדרות מסודרות וסופיות
    candidates = get_ordered_series(all_series)
    if not candidates:
        raise ValueError("No ordered finite series available")

    series_entry = rng.choice(candidates)
    puzzle_func = rng.choice(puzzle_registry)

    # מנסה להעביר גם את כל הסדרות לפונקציה במידה והיא דורשת זאת
    try:
        return puzzle_func(series_entry, rng, all_series=all_series)
    except TypeError:
        try:
            return puzzle_func(series_entry, rng)
        except TypeError:
            # פונקציות שמקבלות רק rng
            return puzzle_func(rng)
