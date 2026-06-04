# Deploy as a web app — free, shareable, password-protected

This guide puts your Polymer Injectivity Calculator online with:
- a **public URL** you can share with anyone,
- a **password** they must enter to use it,
- **zero hosting cost** (Streamlit Community Cloud free tier).

Total time: ~15 minutes.

---

## How it works (the short version)

Your physics engine (`engine.py`) is pure Python — it runs on a server unchanged. `streamlit_app.py` wraps it in a web interface. Streamlit Community Cloud pulls your code from GitHub, runs it, and gives you a URL like `https://eppok-polymer.streamlit.app`. Anyone with the link and password can use it from a browser — no install on their side.

Nothing is stored on the server. Each user uploads their Excel for that session only; it's discarded when they close the tab.

---

## What you need

1. The `eppok-eor-tools` repository pushed to GitHub (see `GITHUB_SETUP.md`).
2. A free [Streamlit Community Cloud](https://share.streamlit.io) account — sign in with your GitHub account.

---

## Step 1 — Test it locally first (optional but recommended)

```bash
cd eppok-eor-tools
pip install -r requirements.txt

# Set a local password
cp .streamlit/secrets.toml.template .streamlit/secrets.toml
# edit .streamlit/secrets.toml and set app_password = "..."

streamlit run streamlit_app.py
```

Your browser opens at `http://localhost:8501`. You'll see the password screen, then the calculator.

---

## Step 2 — Push to GitHub

If you haven't already pushed the repo:

```bash
git add .
git commit -m "feat: add Streamlit web app with password gate and Excel upload"
git push
```

**Important:** `.streamlit/secrets.toml` is gitignored — your real password will NOT be uploaded. That's intentional. You'll set the password directly on Streamlit Cloud in Step 4.

---

## Step 3 — Create the app on Streamlit Cloud

1. Go to [https://share.streamlit.io](https://share.streamlit.io) and sign in with GitHub.
2. Click **Create app** → **Deploy a public app from GitHub**.
3. Fill in:
   - **Repository:** `eppokcompany/eor-tools`
   - **Branch:** `main`
   - **Main file path:** `streamlit_app.py`
   - **App URL:** choose a custom subdomain, e.g. `eppok-polymer` → gives `https://eppok-polymer.streamlit.app`
4. Click **Deploy**. First build takes 2–4 minutes (it installs `requirements.txt`).

> Note on "public app": the app URL is reachable by anyone, but the password gate inside the app blocks use until they enter the password. This is what gives you "share a link, they need the password."

---

## Step 4 — Set the password (Secrets)

1. Once deployed, open your app's menu (⋮ top-right) → **Settings** → **Secrets**.
2. Paste exactly this (with your real password):

   ```toml
   app_password = "your-strong-password-here"
   ```

3. Click **Save**. The app restarts automatically (~30 s).

Now anyone visiting the URL sees the password screen first.

---

## Step 5 — Share

Send your client/colleague two things:
- The URL: `https://eppok-polymer.streamlit.app`
- The password (share it separately — e.g. by phone or a different channel than the link).

They open the link, enter the password, and use the full calculator: presets, all parameters, Excel rate-schedule upload, plots, CSV/PNG download.

---

## Updating the app later

Any time you `git push` to `main`, Streamlit Cloud automatically rebuilds and redeploys within a minute or two. No manual step needed.

To change the password: Settings → Secrets → edit → Save.

---

## Free tier — what to expect

| | Free tier behaviour |
|---|---|
| Cost | $0 |
| Custom URL | Yes (`*.streamlit.app`) |
| Sleeping | App sleeps after ~12h of no traffic; wakes in ~30 s on next visit |
| Resources | 1 GB RAM — plenty for this analytical tool |
| Private repos | Supported (you sign in with GitHub) |

The only visible tradeoff is the wake-from-sleep delay. For a tool you share occasionally, it's a non-issue.

---

## Alternatives (if you outgrow the free tier)

- **Hugging Face Spaces** — also free, similar workflow, can be set fully private.
- **Render / Railway** — free tiers with always-on options behind a paywall; more control.
- **Pure browser (Pyodide + GitHub Pages)** — the engine is pure NumPy, so it *could* run entirely in the browser with no server at all (truly free forever, no sleeping). This is a bigger rebuild — worth it only if the sleep delay becomes annoying or you want zero backend.

For now, Streamlit Community Cloud is the simplest path that meets every requirement: free, nice, link-shareable, password-protected.

---

## Security note (read once)

The shared-password gate is appropriate here because the app stores no data — users upload an Excel per session and it's discarded. The password keeps casual visitors out; it is not designed to protect highly sensitive data sitting on the server (there isn't any). If a client ever needs stricter access control (per-user logins, audit trail), that's a different deployment — ask and we'll scope it.
