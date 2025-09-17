# Firefly Code ✨  

**Firefly** is a publicly available synchrotron radiation code designed for modeling afterglows of high-energy astrophysical transients.  

The synchrotron physics implemented is based on the standard afterglow theory of  
[Sari, Piran & Narayan (1998)](https://doi.org/10.1086/311269).  
Firefly takes a **hydrodynamically evolved astrophysical outflow** and post-processes it to generate simulated observations.  
The version available here includes pre-existing outputs for a **short gamma-ray burst (sGRB) jet** generated using *The Jet Code*.  

Firefly can simulate both **on-axis and off-axis emission** from structured astrophysical outflows, including relativistic jets, across the entire electromagnetic spectrum.  
With Firefly, you can generate:  

- 📈 Multi-wavelength light curves  
- 🌌 Spectra across radio → X-ray  
- 🛰️ 2D sky projections of transients  

over the full dynamical range allowed by the hydrodynamical inputs.  

---

## 🔭 Applications  

Firefly is a versatile tool for interpreting multi-wavelength observations of non-thermal emission from:  

- Gamma-ray bursts (GRBs) and their afterglows  
- Jet–interstellar medium (ISM) interactions  
- Tidal Disruption Events (TDEs)  
- Other transients powered by synchrotron radiation  

---

## ⚖️ Validation & Benchmarking  

Firefly has been standardized and validated against other afterglow codes, including:  

- [BoxFit](https://doi.org/10.1088/0004-637X/751/2/155) (van Eerten & MacFadyen 2012)  
- [pyafterglow / Afterglowpy](https://doi.org/10.3847/1538-4357/ab91b8) (Ryan et al. 2020)  
- [JetSimPy](https://ui.adsabs.harvard.edu/abs/2023ApJ...945...60N) (Nava et al. 2023)  

---

## 🧩 Scientific Use  

Firefly has already been used extensively to study multimessenger and transient events such as:  

- **GRB 170817A** (binary neutron star merger)  
- **PTF10tqv**  
- **EP241021A**  
- Other transients detected by the **Palomar Transient Factory** and **Einstein Probe**  

---

## 📜 Citation  

If you use **Firefly** in your research, please cite:  

> **Dastidar, R. G.** et al. (2025), *Firefly: A Synchrotron Radiation Code for High-Energy Astrophysics*, in preparation.  

### BibTeX  
```bibtex
@misc{dastidar2025firefly,
  author       = {Ranadeep G. Dastidar and collaborators},
  title        = {Firefly: A Synchrotron Radiation Code for High-Energy Astrophysics},
  year         = {2025},
  note         = {in preparation},
  url          = {https://github.com/YOUR-USERNAME/Firefly}
}
