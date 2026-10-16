"""Build a consolidated Excel workbook from the per-tree table: biomass (AGB),
carbon and CO2-equivalent added per tree, with a wood-density sensitivity and a
Colombia-specific cross-check.

Equations, constants and citations come from the primary-source research in
docs/research/biomass-allometry/README.md.
"""
import csv, numpy as np
import openpyxl
from openpyxl.styles import Font, PatternFill

O = "/tmp/opencode/dbh/out"
SRC = f"{O}/tabla_arboles.csv"
OUT = "/home/usuario/workspaces/NubePuntos/prototype/stem-isolation-dbh/results/biomasa_carbono.xlsx"

# ---------------------------------------------------------------- CONFIG
RHO = 0.60          # generic tropical Tier-1 default (species unknown)
RHO_ALT = 0.63      # tropical South America mean (Zanne et al. 2009)
CARBON_FRACTION = 0.47    # IPCC 2006 GL, Vol.4 Ch.4 Table 4.3
CO2_PER_C = 44.0 / 12.0   # = 3.6667


def chave2014(d_cm, h_m, rho):
    """Chave et al. (2014) Eq.4, pantropical DBH+height. kg."""
    return 0.0673 * (rho * d_cm ** 2 * h_m) ** 0.976


def alvarez2012_typeI(d_cm, h_m, rho):
    """Alvarez et al. (2012) Colombia, Type I (D+H+rho). kg.
    ln(AGB) = -2.277 + 2.017 ln D + 0.715 ln H + 0.718 ln rho."""
    return np.exp(-2.277 + 2.017 * np.log(d_cm) + 0.715 * np.log(h_m) + 0.718 * np.log(rho))


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    rows = list(csv.DictReader(open(SRC)))
    for r in rows:
        d = f(r["DBH_cm"]); h = f(r["altura_m"])
        ok = r["estado"] == "measured" and d and h and d > 0 and h > 0
        if ok:
            agb = float(chave2014(d, h, RHO))
            agb_alt = float(chave2014(d, h, RHO_ALT))
            agb_col = float(alvarez2012_typeI(d, h, RHO))
        else:
            agb = agb_alt = agb_col = None
        r["biomasa_aerea_kg"] = round(agb, 1) if agb is not None else ""
        r["carbono_kg"] = round(agb * CARBON_FRACTION, 1) if agb is not None else ""
        r["CO2eq_kg"] = round(agb * CARBON_FRACTION * CO2_PER_C, 1) if agb is not None else ""
        r["biomasa_rho063_kg"] = round(agb_alt, 1) if agb_alt is not None else ""
        r["biomasa_alvarez2012_kg"] = round(agb_col, 1) if agb_col is not None else ""

    wb = openpyxl.Workbook()
    hf = Font(bold=True, color="FFFFFF"); hfill = PatternFill("solid", fgColor="2F5496")

    # ---- Arboles
    ws = wb.active; ws.title = "Arboles"
    cols = ["id_arbol", "latitud", "longitud", "DBH_cm", "altura_m", "confianza",
            "estado", "arco", "seed", "inclinacion_deg", "contrafuerte", "n_pts",
            "pts_banda", "biomasa_aerea_kg", "carbono_kg", "CO2eq_kg",
            "biomasa_rho063_kg", "biomasa_alvarez2012_kg", "motivo"]
    ws.append(cols)
    for c in ws[1]:
        c.font = hf; c.fill = hfill
    for r in rows:
        ws.append([r.get(k, "") for k in cols])
    ws.freeze_panes = "A2"
    for i, w in enumerate([9, 11, 12, 8, 9, 10, 9, 7, 10, 13, 11, 8, 9, 16, 11, 11, 16, 19, 16], 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = w

    # ---- Resumen
    ws2 = wb.create_sheet("Resumen")
    meas = [r for r in rows if r["estado"] == "measured"]
    def tot(rs, k): return round(sum(float(r[k]) for r in rs if r[k] not in ("", None)), 1)
    def row(label, *vals): ws2.append([label, *vals])
    ws2.append(["Resumen de biomasa y carbono (nube: segmentado_crop.laz, 72.16M pts)"])
    ws2["A1"].font = Font(bold=True, size=13)
    row(""); row("Arboles de dosel detectados", len(rows))
    row("Arboles medidos", len(meas)); row("Arboles rechazados", len(rows) - len(meas))
    row("")
    ws2.append(["Ecuacion primaria", "Chave et al. (2014) Eq.4, DBH+altura"]); 
    ws2.append(["rho primario (g/cm3)", RHO]); ws2.append(["rho alternativo (g/cm3)", RHO_ALT])
    ws2.append(["Contraste Colombia", "Alvarez et al. (2012) Type I"]); 
    ws2.append(["Fraccion de carbono", CARBON_FRACTION]); ws2.append(["CO2/C", round(CO2_PER_C, 4)])
    row("")
    ws2.append(["", "biomasa_kg", "carbono_kg", "CO2eq_kg", "CO2eq_t"])
    for c in ws2[ws2.max_row]: c.font = Font(bold=True)
    ws2.append(["Chave 2014 (rho=0.60)", tot(meas, "biomasa_aerea_kg"), tot(meas, "carbono_kg"),
                tot(meas, "CO2eq_kg"), round(tot(meas, "CO2eq_kg") / 1000, 2)])
    ws2.append(["Chave 2014 (rho=0.63)", tot(meas, "biomasa_rho063_kg"),
                round(tot(meas, "biomasa_rho063_kg") * CARBON_FRACTION, 1),
                round(tot(meas, "biomasa_rho063_kg") * CARBON_FRACTION * CO2_PER_C, 1),
                round(tot(meas, "biomasa_rho063_kg") * CARBON_FRACTION * CO2_PER_C / 1000, 2)])
    ws2.append(["Alvarez 2012 (Colombia, rho=0.60)", tot(meas, "biomasa_alvarez2012_kg"),
                round(tot(meas, "biomasa_alvarez2012_kg") * CARBON_FRACTION, 1),
                round(tot(meas, "biomasa_alvarez2012_kg") * CARBON_FRACTION * CO2_PER_C, 1),
                round(tot(meas, "biomasa_alvarez2012_kg") * CARBON_FRACTION * CO2_PER_C / 1000, 2)])
    row("")
    ws2.append(["Por confianza", "n", "biomasa_kg", "CO2eq_kg", "CO2eq_t"])
    for c in ws2[ws2.max_row]: c.font = Font(bold=True)
    for conf in ("high", "medium", "low"):
        rs = [r for r in meas if r["confianza"] == conf]
        ws2.append([conf, len(rs), tot(rs, "biomasa_aerea_kg"), tot(rs, "CO2eq_kg"),
                    round(tot(rs, "CO2eq_kg") / 1000, 3)])
    for i, w in enumerate([36, 14, 14, 14, 12], 1):
        ws2.column_dimensions[openpyxl.utils.get_column_letter(i)].width = w

    # ---- Metodologia
    ws3 = wb.create_sheet("Metodologia")
    txt = [
        "Biomasa aerea, carbono y CO2 equivalente",
        "",
        "Ecuacion primaria (Chave et al. 2014, Eq.4, pantropical, DBH + altura):",
        "  AGB [kg] = 0.0673 * (rho * DBH[cm]^2 * H[m])^0.976",
        "  Fuente: Chave J. et al. (2014), Global Change Biology 20(10):3177-3190, doi:10.1111/gcb.12629.",
        "  Ajustada en 4,004 arboles tropicales talados; CV a nivel de arbol ~56.5%, sesgo ~5.3%.",
        "",
        "rho (densidad de madera especifica, g/cm3):",
        f"  {RHO} = valor por defecto tropical generico (IPCC/FAO Tier-1), especie desconocida.",
        f"  {RHO_ALT} = media de Sudamerica tropical (Zanne et al. 2009, Ecol. Lett. 12:351-366).",
        "  AGB ~ rho^0.976 -> el error de rho se propaga casi 1:1; incertidumbre de rho 10-40%.",
        "",
        "Contraste especifico de Colombia (Alvarez et al. 2012, Type I, D+H+rho):",
        "  ln(AGB) = -2.277 + 2.017 ln(D) + 0.715 ln(H) + 0.718 ln(rho)",
        "  Fuente: Alvarez E. et al. (2012), Forest Ecology and Management 267:297-308.",
        "",
        "Carbono y CO2:",
        f"  C [kg]    = AGB * {CARBON_FRACTION}   (fraccion de carbono por defecto IPCC 2006, Vol.4 Ch.4 Tabla 4.3)",
        f"  CO2eq [kg]= C * 44/12 = AGB * {round(CARBON_FRACTION*CO2_PER_C,4)}   (relacion molar IPCC/UNFCCC)",
        "",
        "Advertencias:",
        "- Sin verdad de campo ni especie: rho es un valor por defecto, no medido.",
        "- Error de ecuacion irreducible (~+-40-60% a nivel de arbol); el error del DBH/altura",
        "  fotogrametricos se propaga (AGB ~ D^1.952: +5% en DBH ~ +10% en AGB).",
        "- Las ecuaciones pantropicales se ajustaron en bosque natural; en agroforesteria los arboles",
        "  de sombra suelen ser mas bajos y anchos para un mismo DBH (posible sesgo).",
        "- Solo se calcula biomasa para arboles 'measured'; los 'refused' quedan vacios.",
        "- Reportar estratificado por confianza, no como un unico numero.",
        "",
        "Documento completo de respaldo (citas, ecuaciones alternativas, unidades):",
        "  docs/research/biomass-allometry/README.md",
    ]
    for line in txt:
        ws3.append([line])
    ws3["A1"].font = Font(bold=True, size=13)
    ws3.column_dimensions["A"].width = 115

    wb.save(OUT)
    print("wrote", OUT)
    print("measured", len(meas),
          "| Chave2014 rho0.60 AGB t", round(tot(meas, "biomasa_aerea_kg") / 1000, 1),
          "CO2eq t", round(tot(meas, "CO2eq_kg") / 1000, 2),
          "| rho0.63 CO2eq t", round(tot(meas, "biomasa_rho063_kg") * CARBON_FRACTION * CO2_PER_C / 1000, 2),
          "| Alvarez CO2eq t", round(tot(meas, "biomasa_alvarez2012_kg") * CARBON_FRACTION * CO2_PER_C / 1000, 2))


if __name__ == "__main__":
    main()
