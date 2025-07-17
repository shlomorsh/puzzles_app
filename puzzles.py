# puzzles.py
# -*- coding: utf-8 -*-
"""אוסף פונקציות חידה + Registry אוטומטי"""

import random
from typing import Dict, Callable, List, Optional

# Registry: כל פונקציית חידה שנרשמת תתווסף לכאן
puzzle_registry: List[Callable] = []


def register(func: Callable) -> Callable:
    """Decorator שמוסיף את הפונקציה למאגר"""
    puzzle_registry.append(func)
    return func


# 'השלם את סימן השאלה' - לוקח סדרה מסודרת,
# בוחר מיקום אות חוקי ומסתיר אות אחת אקראית.
@register
def letter_position_puzzle(series: Dict, rng: random.Random) -> Dict:
    if not series.get("ordered"):
        raise ValueError("Series must be ordered")

    items = series["items"]
    min_len = min(len(w) for w in items)
    position = rng.randrange(min_len)

    letters = [w[position] for w in items]
    missing_index = rng.randrange(len(letters))
    answer = letters[missing_index]
    letters[missing_index] = "?"

    question = ", ".join(letters)
    direction = "התחלה" if position == 0 else "מיקום קבוע בכל מילה"
    hint = f"האות ה-{position + 1} מ{direction}"

    return {
        "title": "השלם את סימן השאלה",
        "question": question,
        "hint": hint,
        "answer": answer,
        "meta": {
            "series_id": series["id"],
            "position": position,
            "missing_index": missing_index,
        },
    }


# 'מצא את ה-?? (שתי אותיות)' - בוחרים slice_mode:
#   - 'first2', 'last2', 'extremes'
# מסתיר זוג אותיות אקראי במערך.
@register
def two_letters_puzzle(series: Dict, rng: random.Random) -> Dict:
    if not series.get("ordered"):
        raise ValueError("Series must be ordered")

    items = series["items"]
    slice_mode = rng.choice(["first2", "last2", "extremes"])

    def pick_letters(w: str) -> str:
        if slice_mode == "first2":
            return w[:2]
        if slice_mode == "last2":
            return w[-2:]
        return w[0] + w[-1]

    excerpts = [pick_letters(w) for w in items]
    missing_index = rng.randrange(len(excerpts))
    answer = excerpts[missing_index]
    excerpts[missing_index] = "??"

    question = ", ".join(excerpts)
    if slice_mode == "first2":
        hint = "שתי האותיות הראשונות בכל מילה"
    elif slice_mode == "last2":
        hint = "שתי האותיות האחרונות בכל מילה"
    else:
        hint = "אות ראשונה + אחרונה בכל מילה"

    return {
        "title": "מצא את ה-?? (שתי אותיות)",
        "question": question,
        "hint": hint,
        "answer": answer,
        "meta": {
            "series_id": series["id"],
            "slice_mode": slice_mode,
            "missing_index": missing_index,
        },
    }


# 'בחר את האופציה המתאימה לרשימה' - בגרסת אנגרמות
# מבנה multiple choice עם 4 אפשרויות (3 הסחות, 1 נכונה).
@register
def multiple_choice_anagram_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    if all_series is None:
        raise ValueError("all_series is required for distractors")

    items = series["items"]
    question_items = rng.sample(items, k=min(4, len(items)))

    def scramble(w: str) -> str:
        letters = list(w)
        rng.shuffle(letters)
        return "".join(letters)

    scrambled_question = [scramble(w) for w in question_items]
    question = ", ".join(scrambled_question)

    remaining = [w for w in items if w not in question_items]
    if not remaining:
        remaining = items.copy()
    correct_word = rng.choice(remaining)
    correct_scrambled = scramble(correct_word)

    pool = []
    for s in all_series:
        if s["id"] != series["id"]:
            pool.extend(s["items"])
    distractors = rng.sample(pool, k=3)
    scrambled_distractors = [scramble(w) for w in distractors]

    options = scrambled_distractors + [correct_scrambled]
    rng.shuffle(options)

    hint = f"המילה הנכונה היא אנגרמה של פריט מ-{series['labels'][0]}"

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "options": options,
        "answer": correct_scrambled,
        "hint": hint,
        "meta": {
            "series_id": series["id"],
            "question_items": question_items,
            "correct_word": correct_word,
        },
    }




# 'השלם את סימן השאלה' – בוחרת סדרה מתאימה עם מספיק מילים מאותו אורך,
# בונה מטריצת תווים (transpose), מסתירה אות אחת אקראית ומחזירה את הרצף.
@register
def matrix_transpose_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    # מוודא שיש רשימת סדרות עבור בחירה חוזרת
    if all_series is None:
        raise ValueError("all_series is required for series fallback")

    # פונקציה למציאת כל הקבוצות של מילים באותו אורך
    def groups_by_length(words):
        lengths = {}
        for w in words:
            lengths.setdefault(len(w), []).append(w)
        return lengths

    # מנסה לספור על הסדרה הנוכחית, אם אינה מספקת – בוחרים סדרה חדשה
    lengths = groups_by_length(series["items"])
    # מצרך: קבוצה של לפחות 3–4 מילים
    valid_lengths = [L for L, grp in lengths.items() if len(grp) >= 3]
    if not valid_lengths:
        # בחר סדרה אחרת עד שתמצא מתאימה
        candidates = [s for s in all_series if s["id"] != series["id"]]
        rng.shuffle(candidates)
        for s in candidates:
            lengths = groups_by_length(s["items"])
            valid_lengths = [L for L, grp in lengths.items() if len(grp) >= 3]
            if valid_lengths:
                series = s
                break
        else:
            raise ValueError("No series with enough words of equal length")

    # בטוח שיש אורך מתאים
    chosen_len = rng.choice(valid_lengths)
    pool = [w for w in series["items"] if len(w) == chosen_len]

    # בוחר 3–4 מילים מהרשימה
    count = min(4, len(pool))
    sel = rng.sample(pool, k=count)

    # בונה מטריצת תווים בתצורת transpose
    # שורה i במטריצה = האות ה-i מכל מילה ב־sel
    letters = []
    for i in range(chosen_len):
        for w in sel:
            letters.append(w[i])

    # מסתיר אות אקראית
    idx = rng.randrange(len(letters))
    answer = letters[idx]
    letters[idx] = "?"

    question = ",".join(letters)
    hint = f"כל המילים בחרו מאורך {chosen_len}; בניית מטריצה (transpose) לפי מיקום תו."

    return {
        "title": "השלם את סימן השאלה",
        "question": question,
        "hint": hint,
        "answer": answer,
        "meta": {
            "series_id": series["id"],
            "selected_words": sel,
            "chosen_length": chosen_len,
            "hidden_index": idx
        },
    }


# 'בחר את האופציה המתאימה לרשימה' – לוקח רשימה Ordered,
# משנה את מיקום הפסיקים בתוך המילים הראשונות, ממשיך עם המילה הבאה כתשובה.
@register
def mispunctuation_word_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    if not series.get("ordered"):
        raise ValueError("Series must be ordered for this puzzle type")
    items = series["items"]
    # בחר 3 מילים רצופות מהסדרה
    start_idx = rng.randrange(0, len(items) - 3)
    seq = items[start_idx:start_idx + 3]
    # הפסקה שגויה: ערבוב מיקומי הפסיק
    def scramble_punct(w: str) -> str:
        # מחלק את המילה לחתיכות אקראיות בין 1 ל3 תווים
        parts = []
        i = 0
        while i < len(w):
            step = rng.randint(1, min(3, len(w) - i))
            parts.append(w[i:i + step])
            i += step
        return ",".join(parts)

    scrambled_seq = [scramble_punct(w) for w in seq]
    # התשובה – המילה הבאה בסדרה
    correct = items[start_idx + 3]
    scrambled_correct = scramble_punct(correct)

    # בחר 3 הסחות - מילים נוספות מהסדרה עם הפסקה שגויה
    distractors = []
    others = [w for i, w in enumerate(items) if i < start_idx or i > start_idx + 3]
    for w in rng.sample(others, k=3):
        distractors.append(scramble_punct(w))

    options = distractors + [scrambled_correct]
    rng.shuffle(options)

    question = ",".join(scrambled_seq)
    hint = "שימרו על סדר ההופעה; ההמשך הוא המילה הרביעית עם אותה 'הפסקה'"

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "options": options,
        "answer": scrambled_correct,
        "hint": hint,
        "meta": {"series_id": series["id"], "start_idx": start_idx},
    }


@register
def country_anagram_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'איזו מדינה זו?' – יוצר אנגרמה (ערבוב אותיות) של שם מדינה מתוך series.json
    ומבקש לזהות איזו מדינה הוטמנה.
    """
    if all_series is None:
        raise ValueError("all_series is required for country selection")

    # מוצא את הסדרה עם id == "countries"
    country_series = next((s for s in all_series if s["id"] == "countries"), None)
    if not country_series:
        raise ValueError("Countries series not found")

    countries = country_series["items"]
    # בוחר מדינה אקראית
    country = rng.choice(countries)

    # מערבב את האותיות עד שלא שווה למקור
    def scramble_word(w: str) -> str:
        letters = list(w)
        scrambled = letters.copy()
        while True:
            rng.shuffle(scrambled)
            if scrambled != letters:
                break
        return "".join(scrambled)

    scrambled = scramble_word(country)
    question = scrambled

    hint = "אנגרמה של שם מדינה"

    return {
        "title": "איזו מדינה זו?",
        "question": question,
        "hint": hint,
        "answer": country,
        "meta": {
            "series_id": country_series["id"],
            "country": country
        },
    }


@register
def letter_shift_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'בחר את האופציה המתאימה לרשימה' – לוקח מספר מילים מתוך series,
    מחליף כל אות באות שלפניה (או אחריה) באלף-בית, ומבקש לנחש איזה מילה שייכת לסדרה.
    """
    if all_series is None:
        raise ValueError("all_series is required for distractors")

    # הגדרת האלף-בית העברי המלא (כולל אותיות סופיות)
    alefbet = list("אבגדהוזחטיכלמנסעפצקרשת")
    final_letters = {"כ": "ך", "מ": "ם", "נ": "ן", "פ": "ף", "צ": "ץ"}
    alefbet += list(final_letters.values())
    shift = rng.choice([-1, 1])  # -1 = קודמת, +1 = עוקבת

    # בוחרים 3 מילים מהסדרה
    items = series["items"]
    q_items = rng.sample(items, k=min(3, len(items)))
    
    # פונקציה להחליף אות
    def shift_word(w: str) -> str:
        out = []
        for ch in w:
            # המרה בין אות סופית לאות רגילה
            ch_norm = ch
            for reg, final in final_letters.items():
                if ch == final:
                    ch_norm = reg
            if ch_norm in alefbet:
                idx = alefbet.index(ch_norm) + shift
                # גלישה בקצוות
                idx = idx % len(alefbet)
                out.append(alefbet[idx])
            else:
                out.append(ch)
        return "".join(out)
    
    # בונים את שאלת הרצף
    scrambled_q = [shift_word(w) for w in q_items]
    question = ",".join(scrambled_q)

    # בוחרים מילה נכונה נוספת מהסדרה, מחוץ ל-q_items
    remaining = [w for w in items if w not in q_items]
    if not remaining:
        remaining = items.copy()
    correct = rng.choice(remaining)
    correct_scr = shift_word(correct)

    # בוחרים 3 הסחות מרשימות אחרות
    pool = [x for s in all_series if s["id"] != series["id"] for x in s["items"]]
    dis = rng.sample(pool, k=3)
    dis_scr = [shift_word(w) for w in dis]

    options = dis_scr + [correct_scr]
    rng.shuffle(options)

    hint_dir = "לפני" if shift == -1 else "אחרי"
    hint = f"האותיות הוחלפו באות שלפניה/אחריה באלף-בית ({hint_dir})"

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "options": options,
        "answer": correct_scr,
        "hint": hint,
        "meta": {
            "series_id": series["id"],
            "shift": shift,
            "q_items": q_items,
            "correct": correct
        },
    }


@register
def mysterious_equation_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'השלם את סימן השאלה' – יוצר רשימת שוויונות שבה צד ימין נגזר מצד שמאל
    לפי כלל אחד מתוך ארבעה (סכומי ספרות, מכפלה+סכום, ריבועים, הפרשים).
    """
    rules = ["sum_pairs", "fixed_sum", "prod_sum", "squares", "abs_diff"]
    rule = rng.choice(rules)

    equations = []
    answers = []

    # מחוללי צד ימין
    def sum_pairs(num):
        a, b, c = [int(d) for d in f"{num:03d}"]
        return int(f"{a+b}{a+c}{b+c}")

    def fixed_sum(num):
        target = 40
        a, b = divmod(num, 10)
        return target - (a + b)

    def prod_sum(num):
        a, b = divmod(num, 10)
        return int(f"{a*b}{a+b}")

    def squares(num):
        a, b = divmod(num, 10)
        return int(f"{a*a}{b*b}")

    def abs_diff(num):
        a, b, c = [int(d) for d in f"{num:03d}"]
        return int(f"{abs(a-b)}{abs(a-c)}")

    func_map = {
        "sum_pairs": sum_pairs,
        "fixed_sum": fixed_sum,
        "prod_sum": prod_sum,
        "squares": squares,
        "abs_diff": abs_diff,
    }

    # בונים 4 שוויונות
    for _ in range(3):
        num = rng.randint(10, 999)
        right = func_map[rule](num)
        equations.append(f"{num} = {right}")

    # השוויון הרביעי עם סימן שאלה
    num_missing = rng.randint(10, 999)
    answer = func_map[rule](num_missing)
    equations.append(f"{num_missing} = ?")

    question = ", ".join(equations)

    return {
        "title": "השלם את סימן השאלה",
        "question": question,
        "hint": "מצא את הכלל הקבוע בין צד שמאל לימין בכל השוויונות",
        "answer": str(answer),
        "meta": {"rule": rule, "missing_left": num_missing},
    }


# numeric_sequence_missing_puzzle (הגרסה התקינה)
@register
def numeric_sequence_missing_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    seq_type = rng.choice(["arithmetic", "geometric", "fibonacci"])
    length = 6
    if seq_type == "arithmetic":
        start = rng.randint(1, 20)
        diff = rng.randint(1, 9)
        seq = [start + i * diff for i in range(length)]
        rule_desc = f"סדרה חשבונית (דלתא {diff})"
    elif seq_type == "geometric":
        start = rng.randint(1, 10)
        ratio = rng.choice([2, 3])
        seq = [start * (ratio ** i) for i in range(length)]
        rule_desc = f"סדרה הנדסית (מכפלה {ratio})"
    else:  # fibonacci
        a, b = rng.randint(1, 10), rng.randint(1, 10)
        seq = [a, b]
        while len(seq) < length:
            seq.append(seq[-1] + seq[-2])
        rule_desc = "פיבונאצ'י"
    hide_idx = rng.randrange(length)
    answer = seq[hide_idx]
    display_seq = [str(x) for x in seq]
    display_seq[hide_idx] = "?"
    return {
        "title": "סדרה חסרה – מה המספר החסר?",
        "question": display_seq,
        "answer": answer,
        "hint": f"שימי לב: זו {rule_desc}.",
        "options": None,
        "meta": {
            "seq_type": seq_type,
            "full_sequence": seq,
            "hidden_index": hide_idx
        }
    }


@register
def quirky_equation_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'בחר את האופציה המתאימה לרשימה' – יוצר 3–4 שוויונות עם כלל חישוב נסתר,
    ומבקש לפתור את השוויון האחרון.
    הכללים האפשריים:
      • איחוד ספרות והכפלה/חילוק
      • חיבור הספרות ואחר־כך כפל קבוע
      • התעלמות מהאופרטור וחלוקה ב־2  (58 -> 29) וכו'.
    """
    rules = ["concat_div2", "digit_sum_x3", "digit_prod_plus"]
    rule = rng.choice(rules)

    def apply_rule(a: int, b: int) -> int:
        if rule == "concat_div2":
            return int(f"{a}{b}") // 2
        if rule == "digit_sum_x3":
            return (a + b) * 3
        if rule == "digit_prod_plus":
            return (a * b) + (a + b)
        return 0  # ברירת מחדל

    equations = []
    for _ in range(3):
        a, b = rng.randint(1, 9), rng.randint(1, 9)
        equations.append(f"{a} ? {b} = {apply_rule(a, b)}")

    # התרגיל החסר
    a, b = rng.randint(1, 9), rng.randint(1, 9)
    answer = apply_rule(a, b)
    question = ", ".join(equations + [f"{a} ? {b} = ?"])

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "answer": str(answer),
        "hint": "אותו כלל הסתתר בכל התרגילים הקודמים",
        "options": [str(answer)] + [str(apply_rule(rng.randint(1, 9), rng.randint(1, 9))) for _ in range(3)],
        "meta": {"rule": rule, "last_a": a, "last_b": b},
    }


@register
def digit_count_range_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'השלם את סימן השאלה' – שואל כמה מופעים של ספרה X יש בטווח A–B.
    """
    start = rng.randint(100, 5000)
    end = start + rng.randint(500, 2000)
    digit = rng.randint(0, 9)

    count = sum(str(n).count(str(digit)) for n in range(start, end + 1))

    question = f"באיזה תדירות מופיעה הספרה {digit} בין {start} ל-{end} (כולל)?"
    return {
        "title": "השלם את סימן השאלה",
        "question": question,
        "hint": "ספר את הופעת הספרה בכל המספרים בטווח",
        "answer": str(count),
        "meta": {"start": start, "end": end, "digit": digit},
    }

@register
def digit_sum_product_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'השלם את השוויון' – שתי ספרות: התוצאה היא (חיבור הספרות)(כפל הספרות).
    לדוגמה: 2 3 → 5 6  |  5 8 → 13 40
    """
    def encode(n: int) -> str:
        a, b = map(int, str(n))
        return f"{a+b}{a*b}"

    # מייצרים 4 משוואות – 3 גלויים ואחד חסר
    nums = [rng.randint(11, 98) for _ in range(4)]          # ללא אפס מוביל
    equations = [f"{n} = {encode(n)}" for n in nums[:3]]

    missing = nums[3]
    answer = encode(missing)
    question = ", ".join(equations + [f'{missing} = ?'])

    return {
        "title": "השלם את השוויון",
        "question": question,
        "hint": "חבר את הספרות ואחר-כך כפל אותן – חבר את התוצאות כמחרוזת אחת",
        "answer": answer,
        "meta": {"rule": "sum_then_product", "missing": missing},
    }


@register
def digit_ops_combo_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'השלם את השוויון – גרסה משוגעת'  
    בוחר אקראית שתי פעולות על הספרות (⊕ , ⊗) ומצמיד את התוצאות.
    הפעולות האפשריות: סכום, מכפלה, הפרש מוחלט, ריבוע-סכום.
    לדוגמה: 4 1  (⊕=סכום, ⊗=הפרש)  → 5 3
    """
    # מגדירים את הפעולות האפשריות
    ops = {
        "sum":       lambda a, b: a + b,
        "product":   lambda a, b: a * b,
        "diff":      lambda a, b: abs(a - b),
        "sum_sq":    lambda a, b: (a + b) ** 2,
    }
    op1_name, op2_name = rng.sample(list(ops.keys()), 2)
    op1, op2 = ops[op1_name], ops[op2_name]

    def encode(n: int) -> str:
        a, b = map(int, str(n))
        return f"{op1(a, b)}{op2(a, b)}"

    # שלוש משוואות לדוגמה ואחת חסרה
    nums = [rng.randint(11, 98) for _ in range(4)]
    equations = [f"{n} = {encode(n)}" for n in nums[:3]]

    missing = nums[3]
    answer = encode(missing)
    question = ", ".join(equations + [f"{missing} = ?"])

    # הסחה קטנה – אפשר להציע גם תשובות שגויות (לא חובה)
    distractors = {encode(rng.randint(11, 98)) for _ in range(6)}
    distractors.discard(answer)
    options = rng.sample(list(distractors), 3) + [answer]
    rng.shuffle(options)

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "options": options,
        "hint": f"שתי פעולות מסתתרות: {op1_name} / {op2_name} (בסדר הזה)",
        "answer": answer,
        "meta": {"op1": op1_name, "op2": op2_name, "missing": missing},
    }

@register
def digit_sum_plus_product_puzzle(series: Dict, rng: random.Random, all_series: Optional[List[Dict]] = None) -> Dict:
    """
    'השלם את השוויון' – שתי ספרות: סכום הספרות + מכפלת הספרות.
    לדוגמה: 5 8 → (5+8) + (5*8) = 13 + 40 = 53
    """
    def encode(n: int) -> int:
        a, b = map(int, str(n))
        return (a + b) + (a * b)

    # שלוש משוואות מלאות ואחת חסרה
    nums = [rng.randint(11, 98) for _ in range(4)]   # נמנעים מאפסים מובילים
    examples = [f"{n} = {encode(n)}" for n in nums[:3]]

    missing = nums[3]
    answer = encode(missing)
    question = ", ".join(examples + [f"{missing} = ?"])

    return {
        "title": "השלם את השוויון",
        "question": question,
        "hint": "חבר את הספרות, כפל אותן, ואז חבר את שתי התוצאות",
        "answer": str(answer),
        "meta": {"rule": "sum_plus_product", "missing_number": missing},
    }


# -------------------------------------------------
# 1. סדרת מספרים חסרה / המשך סדרה
# -------------------------------------------------
@register
def numeric_sequence_missing_puzzle(series: Dict, rng: random.Random,
                                    all_series: Optional[List[Dict]] = None) -> Dict:
    """
    יוצרת רצף בן 6 איברים עם איבר חסר אחד (?). הסוג נבחר אקראית.
    """
    seq_type = rng.choice(["arithmetic", "geometric", "fibonacci"])
    length = 6

    if seq_type == "arithmetic":
        start = rng.randint(1, 20)
        diff = rng.randint(1, 9)
        seq = [start + i * diff for i in range(length)]
        rule_desc = f"סדרה חשבונית (דלתא {diff})"
    elif seq_type == "geometric":
        start = rng.randint(1, 10)
        ratio = rng.choice([2, 3])
        seq = [start * (ratio ** i) for i in range(length)]
        rule_desc = f"סדרה הנדסית (מכפלה {ratio})"
    else:  # fibonacci
        a, b = rng.randint(1, 10), rng.randint(1, 10)
        seq = [a, b]
        while len(seq) < length:
            seq.append(seq[-1] + seq[-2])
        rule_desc = "פיבונאצ'י"

    hide_idx = rng.randrange(length)
    answer = seq[hide_idx]
    display_seq = seq.copy()
    display_seq[hide_idx] = None  # במקום "?"

    return {
        "title": "סדרה חסרה – מה המספר החסר?",
        "question": display_seq,
        "answer": answer,
        "hint": f"שימי לב: זו {rule_desc}.",
        "options": None,
        "meta": {
            "seq_type": seq_type,
            "full_sequence": seq,
            "hidden_index": hide_idx
        }
    }


# -------------------------------------------------
# 2. “משוואות” עם פעולה נסתרת (X ? Y = Z)
# -------------------------------------------------
@register
def hidden_operation_equation_puzzle(series: Dict, rng: random.Random,
                                     all_series: Optional[List[Dict]] = None) -> Dict:
    """
    יוצר 3 “משוואות” עם כלל נסתר, ומבקש להשלים זוג אחרון.
    כלל אפשרי: ax + by + c, כאשר a,b,c קבועים אקראיים קטנים.
    """
    a, b = rng.randint(1, 5), rng.randint(1, 5)
    c = rng.randint(-10, 10)

    def op(x, y):
        return a * x + b * y + c

    # מכינים שלושה זוגות אימון
    equations = []
    used_pairs = set()
    while len(equations) < 3:
        x, y = rng.randint(1, 20), rng.randint(1, 20)
        if (x, y) in used_pairs:
            continue
        used_pairs.add((x, y))
        equations.append((x, y, op(x, y)))

    # זוג הבדיקה
    while True:
        test_x, test_y = rng.randint(1, 20), rng.randint(1, 20)
        if (test_x, test_y) not in used_pairs:
            break
    answer = op(test_x, test_y)

    # טקסט השאלה
    lines = [f"{x} ? {y} = {z}" for x, y, z in equations]
    lines.append(f"{test_x} ? {test_y} = ?")

    return {
        "title": "גלה את הכלל הנסתר",
        "question": "\n".join(lines),
        "answer": answer,
        "hint": "זו פעולה חשבונית ליניארית על x ו-y (לא כפל פשוט).",
        "options": None,
        "meta": {
            "a": a, "b": b, "c": c,
            "training_equations": equations,
            "test_pair": (test_x, test_y)
        }
    }


# -------------------------------------------------
# 3. חידה מילולית־מספרית
# -------------------------------------------------
@register
def word_length_number_puzzle(series: Dict, rng: random.Random,
                              all_series: Optional[List[Dict]] = None) -> Dict:
    """
    נותן רשימת מילים עבריות; המספר הוא מספר האותיות במילה.
    המשתמש צריך למצוא את המספר המתאים למילה האחרונה.
    """
    words = rng.sample([
        "ירושלים", "תל אביב", "חיפה", "באר שבע",
        "אשקלון", "נתניה", "הרצליה", "נהריה",
        "מודיעין", "רחובות", "אילת", "טבריה"
    ], 4)

    pairs = [(w, len(w.replace(" ", ""))) for w in words]  # בלי רווחים
    answer_word, answer_num = pairs[-1]

    question_lines = [f"{w} → {n}" for w, n in pairs[:-1]]
    question_lines.append(f"{answer_word} → ?")

    # מייצר 3 דיסטראקטורים קרובים
    offsets = [-2, -1, 1, 2]
    rng.shuffle(offsets)
    options = [answer_num] + [answer_num + off for off in offsets[:3]]
    rng.shuffle(options)

    return {
        "title": "כמה אותיות יש במילה?",
        "question": "\n".join(question_lines),
        "answer": answer_num,
        "hint": "התוצאה היא פשוט אורך המילה (ללא רווחים).",
        "options": options,
        "meta": {"pairs": pairs}
    }


# -------------------------------------------------
# 4. ספירת תכונה בספרות – כמה ספרות זוגיות?
# -------------------------------------------------
@register
def even_digit_count_puzzle(series: Dict, rng: random.Random,
                            all_series: Optional[List[Dict]] = None) -> Dict:
    nums = []
    while len(nums) < 4:
        n = rng.randint(100, 999)
        if n not in nums:
            nums.append(n)

    def even_count(x):
        return sum(1 for d in str(x) if int(d) % 2 == 0)

    pairs = [(n, even_count(n)) for n in nums]
    answer_num, answer_val = pairs[-1]

    question_lines = [f"{n} → {v}" for n, v in pairs[:-1]]
    question_lines.append(f"{answer_num} → ?")

    options = [answer_val] + rng.sample(
        [k for k in range(4) if k != answer_val], 3
    )
    rng.shuffle(options)

    return {
        "title": "כמה ספרות זוגיות?",
        "question": "\n".join(question_lines),
        "answer": answer_val,
        "hint": "ספר (0,2,4,6,8).",
        "options": options,
        "meta": {"pairs": pairs}
    }


# -------------------------------------------------
# 5. סכום ספרות
# -------------------------------------------------
@register
def digit_sum_puzzle(series: Dict, rng: random.Random,
                     all_series: Optional[List[Dict]] = None) -> Dict:
    nums = rng.sample(range(200, 999), 4)

    def digit_sum(x):
        return sum(int(d) for d in str(x))

    pairs = [(n, digit_sum(n)) for n in nums]
    answer_num, answer_val = pairs[-1]

    lines = [f"{n} → {v}" for n, v in pairs[:-1]]
    lines.append(f"{answer_num} → ?")

    options = [answer_val] + rng.sample(
        [digit_sum(rng.randint(100, 999)) for _ in range(20)
         if digit_sum(_) != answer_val], 3
    )
    rng.shuffle(options)

    return {
        "title": "מה סכום הספרות?",
        "question": "\n".join(lines),
        "answer": answer_val,
        "hint": "פשוט חיבור כל הספרות.",
        "options": options,
        "meta": {"pairs": pairs}
    }


# -------------------------------------------------
# 6. מכפלת ספרות
# -------------------------------------------------
@register
def digit_product_puzzle(series: Dict, rng: random.Random,
                         all_series: Optional[List[Dict]] = None) -> Dict:
    nums = rng.sample(range(112, 987), 4)

    def digit_prod(x):
        prod = 1
        for d in str(x):
            prod *= int(d)
        return prod

    pairs = [(n, digit_prod(n)) for n in nums]
    answer_num, answer_val = pairs[-1]

    lines = [f"{n} → {v}" for n, v in pairs[:-1]]
    lines.append(f"{answer_num} → ?")

    distr = []
    while len(distr) < 3:
        fake = digit_prod(rng.randint(111, 999))
        if fake != answer_val:
            distr.append(fake)
    options = [answer_val] + distr
    rng.shuffle(options)

    return {
        "title": "מכפלת הספרות",
        "question": "\n".join(lines),
        "answer": answer_val,
        "hint": "כפול-כפול-כפול…",
        "options": options,
        "meta": {"pairs": pairs}
    }


# -------------------------------------------------
# 7. שילוב סכום + מכפלה בשרשור
#    דוגמה: 23 → 65  (2+3=5, 2x3=6  → 65)
# -------------------------------------------------
@register
def sum_prod_concat_puzzle(series: Dict, rng: random.Random,
                           all_series: Optional[List[Dict]] = None) -> Dict:
    nums = rng.sample(range(12, 98), 4)  # דו-ספרתי קליל

    def transform(n):
        a, b = divmod(n, 10)
        s = a + b
        p = a * b
        return int(f"{s}{p}")

    pairs = [(n, transform(n)) for n in nums]
    answer_num, answer_val = pairs[-1]

    lines = [f"{n} → {v}" for n, v in pairs[:-1]]
    lines.append(f"{answer_num} → ?")

    options = [answer_val] + rng.sample(
        [transform(rng.randint(12, 98)) for _ in range(30)
         if _ != answer_num and transform(_) != answer_val], 3
    )
    rng.shuffle(options)

    return {
        "title": "חבר ואז כפול – בשרשור",
        "question": "\n".join(lines),
        "answer": answer_val,
        "hint": "חיבור הספרות ואז המכפלה, צמודים אחד לשני.",
        "options": options,
        "meta": {"pairs": pairs}
    }


# -------------------------------------------------
# 8. שרשור ריבועי ספרות
#    דוגמה: 45 → 1625  (4²=16, 5²=25)
# -------------------------------------------------
@register
def square_concat_puzzle(series: Dict, rng: random.Random,
                         all_series: Optional[List[Dict]] = None) -> Dict:
    nums = rng.sample(range(12, 98), 4)

    def transform(n):
        return int("".join(str(int(d) ** 2) for d in str(n)))

    pairs = [(n, transform(n)) for n in nums]
    answer_num, answer_val = pairs[-1]

    lines = [f"{n} → {v}" for n, v in pairs[:-1]]
    lines.append(f"{answer_num} → ?")

    fake_vals = set()
    while len(fake_vals) < 3:
        x = transform(rng.randint(12, 98))
        if x != answer_val:
            fake_vals.add(x)
    options = [answer_val] + list(fake_vals)
    rng.shuffle(options)

    return {
        "title": "שרשור ריבועי ספרות",
        "question": "\n".join(lines),
        "answer": answer_val,
        "hint": "מרובע כל ספרה ומצמידים.",
        "options": options,
        "meta": {"pairs": pairs}
    }


# -------------------------------------------------
# 9. סדרה עם חוק מתחלף:  (-K) , אחר-כך (x-M) , וחוזר חלילה
#    דוגמה מהמאמר: 1,-4,8,3,-6,-11,22,17,...
# -------------------------------------------------
@register
def alternating_rule_sequence_puzzle(series: Dict, rng: random.Random,
                                     all_series: Optional[List[Dict]] = None) -> Dict:
    """
    יוצר סדרה באורך 8-10 שבה החוק מתחלף:
      צעד אי-זוגי: מחסירים K
      צעד זוגי:   מכפילים ב-(‑M)
    """
    length = rng.randint(8, 10)
    first = rng.randint(1, 9)
    K = rng.randint(3, 8)          # כמה מחסירים
    M = rng.choice([2, 3])         # בכמה מכפילים (בסימן שלילי)

    seq = [first]
    sub = True
    while len(seq) < length:
        prev = seq[-1]
        if sub:
            seq.append(prev - K)
        else:
            seq.append(prev * -M)
        sub = not sub  # החלפת החוק

    hide_idx = rng.randrange(3, length - 1)  # מסתירים איבר פנימי
    answer = seq[hide_idx]
    display = [str(x) for x in seq]
    display[hide_idx] = "?"

    return {
        "title": "סדרה עם חוק מתחלף – מצא את המספר החסר",
        "question": display,
        "answer": answer,
        "hint": f"חוק אחד: -{K}, החוק השני: x{-M} – מתחלפים.",
        "options": None,
        "meta": {"K": K, "M": M, "seq": seq, "hidden": hide_idx}
    }


# -------------------------------------------------
# 10. טבלת קלט/פלט ליניארית (ax + b) – להשלים ערכים חסרים
#    מבוסס על הדוגמאות 12→10, 62→35, ...
# -------------------------------------------------
@register
def linear_io_table_puzzle(series: Dict, rng: random.Random,
                           all_series: Optional[List[Dict]] = None) -> Dict:
    """
    יוצר כלל y = ax + b עם a,b שלמים קטנים.
    מציג 6-7 זוגות, חלקם חסרים, ומבקש למלא.
    """
    a = rng.choice([2, 3, 4, -2, -3])
    b = rng.randint(-10, 10)

    n_rows = 6
    inputs = rng.sample(range(10, 90), n_rows)
    rows = []
    for x in inputs:
        rows.append([x, a * x + b])

    # מחביאים 2-3 יציאות ואולי כניסה אחת
    missing_out = rng.sample(range(n_rows), 3)
    missing_in  = rng.choice(range(n_rows))

    display_rows = []
    for idx, (x, y) in enumerate(rows):
        disp_x = "?" if idx == missing_in and rng.random() < 0.4 else x
        disp_y = "?" if idx in missing_out else y
        display_rows.append(f"{disp_x} → {disp_y}")

    # האתגר: למצוא y במקום הראשון ברשימת missing_out
    target_idx = missing_out[0]
    answer = rows[target_idx][1]

    return {
        "title": "טבלת קלט-פלט – גלה את הכלל והשלים",
        "question": "\n".join(display_rows),
        "answer": answer,
        "hint": "הכלל ליניארי: y = a·x + b.",
        "options": None,
        "meta": {"a": a, "b": b, "rows": rows,
                 "missing_out": missing_out, "missing_in": missing_in}
    }


# -------------------------------------------------
# 11. סדרה ריבועית:  n² ± c   (נולדה מה-n²-1 שבמאמר)
# -------------------------------------------------
@register
def quadratic_sequence_puzzle(series: Dict, rng: random.Random,
                              all_series: Optional[List[Dict]] = None) -> Dict:
    """
    יוצר סדרה באורך 7 לפי הנוסחה  n² ± c  (c קבוע קטן).
    מסתיר איבר אמצעי.
    """
    c = rng.choice([-3, -1, 1, 2, 5])
    length = 7
    seq = [(n ** 2) + c for n in range(length)]
    hide_idx = rng.randrange(2, 5)
    answer = seq[hide_idx]
    display = [str(x) for x in seq]
    display[hide_idx] = "?"
    sign = "+" if c > 0 else "-"
    abs_c = abs(c)
    return {
        "title": "סדרה ריבועית – מה המספר החסר?",
        "question": display,
        "answer": answer,
        "hint": f"האיבר ה-n הוא n² {sign} {abs_c}.",
        "options": None,
        "meta": {"c": c, "seq": seq, "hidden": hide_idx}
    }


# -------------------------------------------------
# 12. כלל דו-שלבי (x a ואז +b) על רצף מספרים לא רציף
#    כמו 2→7, 3→10, 5→16  (y = 3x+1)
# -------------------------------------------------
@register
def two_step_linear_sequence_puzzle(series: Dict, rng: random.Random,
                                    all_series: Optional[List[Dict]] = None) -> Dict:
    """
    בוחר a,b (a≠0), מייצר 5-6 זוגות x→y לא רציפים.
    משתמש בשניים-שלושה זוגות להדגמת הכלל, מסתיר את היתר.
    """
    a = rng.choice([2, 3, 4, 5])
    b = rng.randint(-5, 6)

    xs = rng.sample(range(1, 15), 6)
    pairs = [(x, a * x + b) for x in xs]

    # מציגים שלושה זוגות מלאים ושניים חסרים
    rng.shuffle(pairs)
    known = pairs[:3]
    hidden = pairs[3:]

    lines = [f"{x} → {y}" for x, y in known]
    for x, _ in hidden:
        lines.append(f"{x} → ?")

    answer = hidden[0][1]

    return {
        "title": "כלל כפול-ואז-חיבור – מצא את ה-y החסר",
        "question": "\n".join(lines),
        "answer": answer,
        "hint": "קודם כופלים במספר קבוע, אחר-כך מוסיפים קבוע נוסף.",
        "options": None,
        "meta": {"a": a, "b": b, "pairs": pairs}
    }

