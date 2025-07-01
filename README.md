# מחולל החידות (Puzzles App)

מערכת אינטראקטיבית ליצירת ופתרון חידות בעברית, מבוססת Streamlit.

## תיאור
- מחולל חידות אוטומטי בעברית
- מגוון סוגי חידות (מילוליות, סדרות, אנגרמות, מספרים ועוד)
- ממשק וובי נוח (Streamlit)
- מבוסס קבצי JSON

## התקנה והרצה
1. ודא שיש לך Python 3.7+
2. התקן את Streamlit:
   ```bash
   pip install streamlit
   ```
3. להרצה:
   ```bash
   streamlit run app.py
   ```
   או להריץ את הקובץ run_app.bat ב-Windows.

## דוגמה
![דוגמה לממשק](screenshot.png)

## קבצים עיקריים
- `app.py` – ממשק המשתמש
- `puzzle_engine.py` – מנוע החידות
- `puzzles.py` – מאגר פונקציות החידות
- `series.json` – בסיס הנתונים של הסדרות
- `PROJECT_DEFINITION.md` – תיעוד ארכיטקטורה
- `DEVELOPMENT_GUIDE.md` – מדריך פיתוח

## תרומה
תרגישו חופשי להציע חידות חדשות, סדרות, או לשפר את הממשק!

## רישיון
MIT License

---

# Puzzles App (מחולל החידות)

An interactive puzzle generator and solver in Hebrew, based on Streamlit.

## Description
- Automatic puzzle generator in Hebrew
- Many puzzle types (verbal, series, anagrams, numbers, etc.)
- Web interface (Streamlit)
- JSON-based data

## Installation & Run
1. Make sure you have Python 3.7+
2. Install Streamlit:
   ```bash
   pip install streamlit
   ```
3. To run:
   ```bash
   streamlit run app.py
   ```
   Or run `run_app.bat` on Windows.

## Main Files
- `app.py` – Web interface
- `puzzle_engine.py` – Puzzle engine
- `puzzles.py` – Puzzle functions
- `series.json` – Data series
- `PROJECT_DEFINITION.md` – Architecture docs
- `DEVELOPMENT_GUIDE.md` – Dev guide

## Contribution
Feel free to suggest new puzzles, series, or UI improvements!

## License
MIT License 