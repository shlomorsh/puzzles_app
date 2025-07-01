# app.py
# -*- coding: utf-8 -*-
import streamlit as st
import random

from puzzle_engine import load_series, generate_random_puzzle

# ----- הגדרות בסיס -----
SEED = None
rng = random.Random(SEED)
all_series = load_series()

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
st.subheader(puzzle["title"])
st.write("**" + puzzle["question"] + "**")

# ----- קלט משתמש (בחירה מרובה או טקסט) -----
if "options" in puzzle:
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
        value=st.session_state.user_answer,
        key="user_input"
    )
    st.session_state.user_answer = user_input

# ----- כפתור רמז -----
if st.button("💡 רמז", key="hint_btn"):
    if not st.session_state.used_hint:
        st.info(puzzle["hint"])
        st.session_state.used_hint = True
    else:
        st.warning("כבר קיבלת רמז 😉")

# ----- בדיקת תשובה -----
if st.button("✅ בדוק תשובה", key="check_btn"):
    answer = st.session_state.user_answer.strip()
    correct = puzzle["answer"]
    if not answer:
        st.warning("כתוב משהו קודם…")
    elif answer == correct:
        st.success("בול! כל הכבוד 🎉")
    else:
        st.error("לא נכון, נסה שוב או בקש רמז.")

# ----- תצוגת פרטים טכניים (אופציונלי) -----
with st.expander("🔍 פרטים טכניים (dev)"):
    st.json(puzzle.get("meta", {}))
