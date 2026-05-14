# How to publish this folder as a GitHub Page

The directory `nc-integers-pages/` is a self-contained GitHub Pages site.
Follow these steps once to publish it; afterwards every `git push` updates
both the site and the certificate archive.

## 1. Create the repository on GitHub

1. Sign in to GitHub.
2. Click *New repository*.
3. Name it `nc-integers-certificates` (or any name; update `index.md` and
   `README.md` accordingly).
4. Visibility: **Public** (required for free GitHub Pages).
5. Do **not** initialize with README / .gitignore / LICENSE — this folder
   already contains them.

## 2. Initialize git locally and push

From inside `nc-integers-pages/`:

```bash
git init -b main
git add .
git commit -m "Initial commit: NC integers supplementary archive"
git remote add origin https://github.com/<your-github-username>/nc-integers-certificates.git
git push -u origin main
```

## 3. Enable GitHub Pages

1. Open the repository on GitHub.
2. *Settings* → *Pages*.
3. Source: **Deploy from a branch**.
4. Branch: **main**, folder: **/ (root)**.
5. Click *Save*.

After a few seconds the site appears at

```
https://<your-github-username>.github.io/nc-integers-certificates/
```

## 4. Update `index.md` and `README.md`

Search and replace `<your-github-username>` in both files with your actual
GitHub handle, then push again:

```bash
git add index.md README.md
git commit -m "Wire in GitHub username"
git push
```

## 5. (Optional) Enable the CI verification workflow

The file `.github/workflows/verify.yml` regenerates the certificate tree on
every push and uploads it as a build artifact. No setup is required beyond
having the file in the repository; the workflow runs automatically.

## 6. Cite

Once the site is live, paste the actual URL into

- the *Computational methods* section of `NC_Integers_v3.tex`, replacing
  the placeholder "supplementary file accompanying this article";
- the `Corradetti-NC-Integers-Certificates-2026` bibtex entry in
  `README.md`.

## Local preview (optional)

To preview the GitHub Page locally with Jekyll:

```bash
gem install bundler jekyll
bundle init
bundle add github-pages jekyll-remote-theme jekyll-relative-links jekyll-optional-front-matter
bundle exec jekyll serve
```

The site appears at <http://127.0.0.1:4000/>.

## Re-running verification locally

```bash
python verify.py --clean
```

Requires Python 3.11. No third-party packages needed.
