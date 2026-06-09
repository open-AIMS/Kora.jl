import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/Kora.jl/dev/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Start Here', link: '/what-can-kora-tell-me' },
{ text: 'Getting Started', link: '/getting-started' },
{ text: 'Model Overview', link: '/model-overview' },
{ text: 'Concepts', collapsed: false, items: [
{ text: 'Decision Support Under Uncertainty', link: '/concepts/decision-support-under-uncertainty' }]
 },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Running Simulations', link: '/tutorials/running-simulations' },
{ text: 'Visualization', link: '/tutorials/visualization' },
{ text: 'Ensemble Analysis', link: '/tutorials/ensemble-analysis' },
{ text: 'Fitting from EcoRRAP', link: '/tutorials/fitting-from-ecorrap' },
{ text: 'Restoration Scenarios', link: '/tutorials/restoration-scenarios' }]
 },
{ text: 'Calibration', collapsed: false, items: [
{ text: 'Model Calibration', link: '/calibration/model-calibration' },
{ text: 'Ensemble Assessment', link: '/calibration/ensemble-assessment' }]
 },
{ text: 'API Reference', collapsed: false, items: [
{ text: 'Reef State', link: '/reference/api-reef-state' },
{ text: 'Simulation', link: '/reference/api-simulation' },
{ text: 'Coral Models', link: '/reference/api-models' },
{ text: 'Model I/O', link: '/reference/api-interface' },
{ text: 'Coral Dynamics', link: '/reference/api-coral-dynamics' },
{ text: 'Metrics', link: '/reference/api-metrics' }]
 },
{ text: 'Background', collapsed: false, items: [
{ text: 'Coral Biology', link: '/concepts/coral-biology' }]
 },
{ text: 'Glossary', link: '/glossary' },
{ text: 'Contributing', link: '/contributing' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/Kora.jl/dev/',// TODO: replace this in makedocs!
  title: 'Kora.jl',
  description: 'Documentation for Kora.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  markdown: {
    codeTransformers: [juliaReplTransformer()],
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('/Kora.jl'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  themeConfig: {
    outline: 'deep',
    
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Start Here', link: '/what-can-kora-tell-me' },
{ text: 'Getting Started', link: '/getting-started' },
{ text: 'Model Overview', link: '/model-overview' },
{ text: 'Concepts', collapsed: false, items: [
{ text: 'Decision Support Under Uncertainty', link: '/concepts/decision-support-under-uncertainty' }]
 },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Running Simulations', link: '/tutorials/running-simulations' },
{ text: 'Visualization', link: '/tutorials/visualization' },
{ text: 'Ensemble Analysis', link: '/tutorials/ensemble-analysis' },
{ text: 'Fitting from EcoRRAP', link: '/tutorials/fitting-from-ecorrap' },
{ text: 'Restoration Scenarios', link: '/tutorials/restoration-scenarios' }]
 },
{ text: 'Calibration', collapsed: false, items: [
{ text: 'Model Calibration', link: '/calibration/model-calibration' },
{ text: 'Ensemble Assessment', link: '/calibration/ensemble-assessment' }]
 },
{ text: 'API Reference', collapsed: false, items: [
{ text: 'Reef State', link: '/reference/api-reef-state' },
{ text: 'Simulation', link: '/reference/api-simulation' },
{ text: 'Coral Models', link: '/reference/api-models' },
{ text: 'Model I/O', link: '/reference/api-interface' },
{ text: 'Coral Dynamics', link: '/reference/api-coral-dynamics' },
{ text: 'Metrics', link: '/reference/api-metrics' }]
 },
{ text: 'Background', collapsed: false, items: [
{ text: 'Coral Biology', link: '/concepts/coral-biology' }]
 },
{ text: 'Glossary', link: '/glossary' },
{ text: 'Contributing', link: '/contributing' }
]
,
    sidebarDrawer: false,
    editLink: { pattern: "https://https://github.com/open-AIMS/Kora.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/open-AIMS/Kora.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
