// MathJax configuration for the manual.
// It loads mhchem, which allows chemical formulas with \ce{...}, for example $\ce{C60}$ or $\ce{C6H6}$.
window.MathJax = {
  loader: {
    load: ["[tex]/mhchem"]
  },
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    packages: { "[+]": ["mhchem"] }
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};
