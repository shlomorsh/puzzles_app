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
        puzzle = puzzle_func(series_entry, rng, all_series=all_series)
    except TypeError:
        try:
            puzzle = puzzle_func(series_entry, rng)
        except TypeError:
            # פונקציות שמקבלות רק rng
            puzzle = puzzle_func(rng)
    # הוספת שם הפונקציה ל-meta
    if "meta" in puzzle and isinstance(puzzle["meta"], dict):
        puzzle["meta"]["func_name"] = puzzle_func.__name__
    else:
        puzzle["meta"] = {"func_name": puzzle_func.__name__}
    return puzzle


def debug_puzzle_coverage(all_series):
    import traceback
    print("\n--- בדיקת זמינות חידות ---")
    for i, func in enumerate(puzzle_registry):
        func_name = func.__name__
        successes = 0
        failures = []
        for s in all_series:
            try:
                # ננסה להריץ עם כל הפרמטרים האפשריים
                try:
                    func(s, random.Random(), all_series=all_series)
                except TypeError:
                    try:
                        func(s, random.Random())
                    except TypeError:
                        func(random.Random())
                successes += 1
            except Exception as e:
                failures.append(f"סדרה: {s.get('id', '?')} | {str(e)}")
        if successes:
            print(f"✅ {func_name}: הצליח עבור {successes} סדרות מתוך {len(all_series)}")
        if failures:
            print(f"❌ {func_name}: נכשל עבור {len(failures)} סדרות:")
            for fail in failures:
                print(f"   - {fail}")
    print("--- סוף בדיקה ---\n")


def list_working_puzzles():
    """
    מדפיסה למסוף את שמות כל החידות שנכשלות (לא עובדות) עם פירוט החריגה, עם סדרה דמה.
    """
    import random
    from puzzles import puzzle_registry
    mock_series = {
        "id": "mock",
        "labels": ["דמה"],
        "items": ["אבא", "אמא", "דוד", "רון", "תמר", "שיר"],
        "ordered": True,
        "finite": True
    }
    mock_all_series = [
        mock_series,
        {"id": "mock2", "labels": ["דמה2"], "items": ["אור", "ים", "טל", "חן", "גל"], "ordered": True, "finite": True},
        {"id": "countries", "labels": ["מדינות"], "items": ["ישראל", "צרפת", "גרמניה", "ספרד", "איטליה"], "ordered": True, "finite": True}
    ]
    print("\n--- חידות שנכשלו (ללא תלות בסדרות) ---")
    any_failed = False
    for func in puzzle_registry:
        try:
            try:
                func(mock_series, random.Random(), all_series=mock_all_series)
            except TypeError:
                try:
                    func(mock_series, random.Random())
                except TypeError:
                    func(random.Random())
        except Exception as e:
            any_failed = True
            print(f"❌ {func.__name__}: {type(e).__name__}: {e}")
    if not any_failed:
        print("(כל החידות נטענות בהצלחה)")
