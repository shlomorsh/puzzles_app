# puzzles.py
# -*- coding: utf-8 -*-
"""אוסף פונקציות חידה + Registry אוטומטי"""

import random
from typing import Dict, Callable, List

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
def multiple_choice_anagram_puzzle(series: Dict, rng: random.Random, all_series: List[Dict] = None) -> Dict:
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


# 'בחר את האופציה המתאימה לרשימה' - מטמיע מילה קצרה מתוך סדרה בתוך מילה ארוכה.
# יוצר 4 אפשרויות, רק אחת מכילה טקסט ממוטמן.
@register
def embedded_word_puzzle(series: Dict, rng: random.Random, all_series: List[Dict] = None) -> Dict:
    if all_series is None:
        raise ValueError("all_series is required for distractors")

    items = series["items"]
    short_word = rng.choice(items)

    mode = rng.choice(["start", "end", "random"])
    def embed(word: str) -> str:
        filler = "אבגדהוזחטיכלמנסעפצקרשת"
        n = rng.randint(3, 6)
        pad = ''.join(rng.choice(list(filler)) for _ in range(n))
        if mode == "start":
            return short_word + pad
        if mode == "end":
            return pad + short_word
        pos = rng.randint(0, len(word))
        return word[:pos] + short_word + word[pos:]

    pool = []
    for s in all_series:
        if s["id"] != series["id"]:
            pool.extend(s["items"])
    distractor_words = rng.sample(pool, k=3)

    options = [embed(short_word)]
    for w in distractor_words:
        options.append(embed(w))
    rng.shuffle(options)

    question = ", ".join(options)
    hint = f"המילה הנכונה כוללת מילת מפתח מתוך {series['labels'][0]}"

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "options": options,
        "answer": next(opt for opt in options if short_word in opt),
        "hint": hint,
        "meta": {
            "series_id": series["id"],
            "mode": mode,
            "short_word": short_word,
        },
    }


# 'השלם את סימן השאלה' – בוחרת סדרה מתאימה עם מספיק מילים מאותו אורך,
# בונה מטריצת תווים (transpose), מסתירה אות אחת אקראית ומחזירה את הרצף.
@register
def matrix_transpose_puzzle(series: Dict, rng: random.Random, all_series: List[Dict] = None) -> Dict:
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
def mispunctuation_word_puzzle(series: Dict, rng: random.Random, all_series: List[Dict] = None) -> Dict:
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
def numeric_sequence_puzzle(rng: random.Random) -> Dict:
    """
    'בחר את האופציה המתאימה לרשימה' – מתחיל ממספר אקראי בן 4 ספרות,
    מפריד לפסקאות אקראיות וממשיך במספר העוקב.
    """
    # התחלה רנדומלית
    start = rng.randint(1000, 9999)
    seq_nums = [start + i for i in range(4)]
    # פונקציית פיצול
    def split_number(n: int) -> str:
        s = str(n)
        parts = []
        i = 0
        while i < len(s):
            step = rng.randint(1, len(s) - i)
            parts.append(s[i:i + step])
            i += step
        return ",".join(parts)

    scrambled_seq = [split_number(n) for n in seq_nums]
    question = ",".join(scrambled_seq)

    correct_num = start + 4
    scrambled_correct = split_number(correct_num)
    # הסחות - מספרים נוספים
    distractors = [split_number(rng.randint(1000, 9999)) for _ in range(3)]
    options = distractors + [scrambled_correct]
    rng.shuffle(options)

    hint = "המשך סדרת המספרים לפי הערך הנוסף עם אותה הפסקת ספרות"

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "options": options,
        "answer": scrambled_correct,
        "hint": hint,
        "meta": {"start_number": start},
    }


@register
def country_anagram_puzzle(series: Dict, rng: random.Random, all_series: List[Dict] = None) -> Dict:
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
def letter_shift_puzzle(series: Dict, rng: random.Random, all_series: List[Dict] = None) -> Dict:
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
def mysterious_equation_puzzle(rng: random.Random) -> Dict:
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


@register
def numeric_sequence_missing_puzzle(rng: random.Random) -> Dict:
    """
    'השלם את סימן השאלה' – יוצר סדרה מספרית עם איבר חסר (שאלה 2).
    בוחר אחד מחמישה סוגי סדרות.
    """
    series_type = rng.choice(
        ["fibo", "arithmetic", "geometric", "digit_product", "cumulative_sum"]
    )

    seq = []
    if series_type == "fibo":
        a, b = rng.randint(1, 20), rng.randint(1, 20)
        seq = [a, b]
        for _ in range(4):
            seq.append(seq[-1] + seq[-2])

    elif series_type == "arithmetic":
        start = rng.randint(1, 50)
        diff = rng.randint(2, 10)
        seq = [start + diff * i for i in range(6)]

    elif series_type == "geometric":
        start = rng.randint(2, 6)
        ratio = rng.randint(2, 4)
        seq = [start * (ratio**i) for i in range(6)]

    elif series_type == "digit_product":
        n = rng.randint(21, 99)
        seq = [n]
        for _ in range(5):
            prod = 1
            for d in str(seq[-1]):
                prod *= int(d)
            seq.append(prod)

    elif series_type == "cumulative_sum":
        vals = [rng.randint(1, 20) for _ in range(6)]
        running = []
        total = 0
        for v in vals:
            total += v
            running.append(total)
        seq = running

    # מסתיר איבר אקראי (לא הראשון)
    idx = rng.randint(1, len(seq) - 2)
    answer = seq[idx]
    seq[idx] = "?"

    question = ", ".join(map(str, seq))
    return {
        "title": "השלם את סימן השאלה",
        "question": question,
        "hint": "גלה את החוקיות של הסדרה ומצא את המספר החסר",
        "answer": str(answer),
        "meta": {"series_type": series_type, "missing_index": idx},
    }


@register
def next_palindrome_number_puzzle(rng: random.Random) -> Dict:
    """
    'בחר את האופציה המתאימה לרשימה' – מוצא את הפלינדרום הבא (או קודם) של מספר בן 3–5 ספרות.
    """
    direction = rng.choice(["next", "prev"])
    num = rng.randint(100, 99999)

    def is_pal(n: int) -> bool:
        s = str(n)
        return s == s[::-1]

    step = 1 if direction == "next" else -1
    candidate = num + step
    while not is_pal(candidate):
        candidate += step

    question = f"{num} – מהו הפלינדרום ה{ 'הבא' if direction=='next' else 'הקודם' }?"
    answer = str(candidate)

    # הסחות – עוד שלושה מספרים אקראיים לא פלינדרומיים
    distractors = []
    while len(distractors) < 3:
        d = rng.randint(num - 500, num + 500)
        if not is_pal(d) and d != candidate:
            distractors.append(str(d))

    options = distractors + [answer]
    rng.shuffle(options)

    return {
        "title": "בחר את האופציה המתאימה לרשימה",
        "question": question,
        "options": options,
        "answer": answer,
        "hint": "מספר שקוראים אותו אותו דבר משמאל לימין ומימין לשמאל",
        "meta": {"direction": direction, "base_number": num},
    }


@register
def quirky_equation_puzzle(rng: random.Random) -> Dict:
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
def digit_count_range_puzzle(rng: random.Random) -> Dict:
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


