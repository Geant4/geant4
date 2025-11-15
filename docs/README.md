# Geant4 Documentation Site

This directory contains the VitePress-based documentation site for Geant4.

## Prerequisites

- Node.js 18+ or 20+
- npm

## Local Development

1. Install dependencies:
   ```bash
   npm install
   ```

2. Start the development server:
   ```bash
   npm run docs:dev
   ```

3. Open your browser to `http://localhost:5173`

## Build

Build the static site:

```bash
npm run docs:build
```

The built site will be in `.vitepress/dist/`

## Preview

Preview the built site locally:

```bash
npm run docs:preview
```

## Deployment

The documentation is automatically deployed to GitHub Pages when changes are pushed to the main/master branch.

The deployment is handled by the GitHub Actions workflow in `.github/workflows/deploy-docs.yml`.

## Site Structure

```
docs/
├── .vitepress/
│   └── config.js          # VitePress configuration
├── reference/             # Reference documentation
│   ├── index.md          # Reference overview
│   ├── build-system.md   # Build system details
│   ├── source-modules.md # Source modules reference
│   └── contributing.md   # Contributing guide
├── index.md              # Home page
├── getting-started.md    # Getting started guide
├── architecture.md       # Architecture overview
└── package.json          # Node.js dependencies
```

## Adding Documentation

1. Create new `.md` files in the appropriate directory
2. Update `.vitepress/config.js` to add navigation/sidebar entries
3. Use standard Markdown with VitePress extensions

## VitePress Features

- Markdown-based content
- Vue components in Markdown
- Built-in search
- Responsive design
- Fast hot module replacement during development

## Learn More

- [VitePress Documentation](https://vitepress.dev/)
- [Geant4 Official Site](https://cern.ch/geant4)
