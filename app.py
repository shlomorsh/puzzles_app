# app.py
# -*- coding: utf-8 -*-
import streamlit as st
import random

from puzzle_engine import load_series, generate_random_puzzle, debug_puzzle_coverage, list_working_puzzles
from puzzles import puzzle_registry

# ----- הגדרות בסיס -----
SEED = None
rng = random.Random(SEED)
all_series = load_series()
debug_puzzle_coverage(all_series)
list_working_puzzles()

# ----- Session state initialization -----
if "current_puzzle" not in st.session_state:
    st.session_state.current_puzzle = None
if "used_hint" not in st.session_state:
    st.session_state.used_hint = False
if "user_answer" not in st.session_state:
    st.session_state.user_answer = ""

def new_puzzle():
    st.session_state.current_puzzle = generate_random_puzzle(rng, all_series)
    st.session_state.used_hint = False
    st.session_state.user_answer = ""

# ----- כותרת -----
st.title("🧩 מחולל החידות – ממשק דפדפן")

# ----- כפתור חידה חדשה -----
if st.button("🎲 חידה חדשה", key="new_puzzle_btn", on_click=new_puzzle):
    pass

# ----- יצירה ראשונית של חידה אם אין אחת -----
if st.session_state.current_puzzle is None:
    new_puzzle()

# ----- וידוא שהמשתנה puzzle קיים -----
puzzle = st.session_state.current_puzzle

# ----- הצגת החידה -----
# הצגת מספר החידה הנוכחית מתוך סך כל החידות
# נמצא את האינדקס של הפונקציה שנבחרה בפועל
puzzle_func_name = None
if puzzle and isinstance(puzzle, dict):
    meta = puzzle.get("meta", {})
    puzzle_func_name = meta.get("func_name")
if puzzle_func_name:
    try:
        current_index = [f.__name__ for f in puzzle_registry].index(puzzle_func_name) + 1
    except ValueError:
        current_index = 0
else:
    current_index = 0
num_puzzles = len(puzzle_registry)
st.markdown(f"<div style='direction:rtl; font-size:18px; color:#888;'>חידה {current_index} מתוך {num_puzzles}</div>", unsafe_allow_html=True)

# אסיר את הצגת כותרת החידה מהממשק
# מחקתי את השורה:
# st.subheader(puzzle["title"] if puzzle and isinstance(puzzle, dict) and "title" in puzzle else "")
st.write("**" + str(puzzle["question"] if puzzle and isinstance(puzzle, dict) and "question" in puzzle else "") + "**")

# ----- קלט משתמש (בחירה מרובה או טקסט) -----
if puzzle and isinstance(puzzle, dict) and "options" in puzzle and puzzle["options"]:
    # Multiple-choice
    user_choice = st.radio(
        "מה התשובה הנכונה?",
        puzzle["options"],
        index=0,
        key="user_choice"
    )
    st.session_state.user_answer = user_choice
else:
    # Free-text input
    user_input = st.text_input(
        "השלם את סימן השאלה:",
        value=st.session_state.user_answer if hasattr(st.session_state, 'user_answer') else "",
        key="user_input"
    )
    st.session_state.user_answer = user_input

# ----- כפתור רמז -----
if st.button("💡 רמז", key="hint_btn"):
    if not st.session_state.used_hint:
        st.info(puzzle["hint"] if puzzle and isinstance(puzzle, dict) and "hint" in puzzle else "")
        st.session_state.used_hint = True
    else:
        st.warning("כבר קיבלת רמז 😉")

# ----- בדיקת תשובה -----
if st.button("✅ בדוק תשובה", key="check_btn"):
    answer = st.session_state.user_answer.strip() if hasattr(st.session_state, 'user_answer') and st.session_state.user_answer else ""
    correct = puzzle["answer"] if puzzle and isinstance(puzzle, dict) and "answer" in puzzle else ""
    if not answer:
        st.warning("כתוב משהו קודם…")
    elif answer == correct:
        st.success("בול! כל הכבוד 🎉")
    else:
        st.error("לא נכון, נסה שוב או בקש רמז.")

# ----- תצוגת פרטים טכניים (אופציונלי) -----
with st.expander("🔍 פרטים טכניים (dev)"):
    st.json(puzzle["meta"] if puzzle and isinstance(puzzle, dict) and "meta" in puzzle else {})
