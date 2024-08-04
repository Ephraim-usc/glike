import os
import pickle
import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mpltern
import glike

font_dirs = ["fonts"]  # The path to the custom font file.
font_files = mpl.font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
  mpl.font_manager.fontManager.addfont(font_file)

plt.rcParams["font.family"] = "Arial"


def adjust_color(color, alpha):
  if type(color) == str:
    color = mcolors.to_rgba(color)
  elif len(color) == 3:
    color = (color[0], color[1], color[2], 1.0)
  color_adj = (color[0], color[1], color[2], color[3]*alpha)
  return color_adj


def boxplot(ax, errors, title, names, colors):
  errors[errors > 3] = 3
  n = errors.shape[0]
  
  ax.axhline(y=0, xmin=0, xmax=1, linestyle = "solid", color = "black", zorder = 1)
  ax.axhline(y=-0.3, xmin=0, xmax=1, linestyle = "dotted", color = "black", zorder = 1)
  ax.axhline(y=0.3, xmin=0, xmax=1, linestyle = "dotted", color = "black", zorder = 1)
  
  for i in range(len(truth)):
    ax.scatter([i+1 - 0.2 + 0.4*np.random.uniform() for _ in range(n)], errors.iloc[:, i], s = 5, facecolors='none', edgecolor = colors[i], linewidths = 0.5, zorder=2)
  
  bp = ax.boxplot(errors.values, showfliers = False, showcaps = False, zorder = 2, widths = 0.75, patch_artist=True)
  for i, color in enumerate(colors):
    bp["boxes"][i].set_edgecolor(colors[i])
    bp["boxes"][i].set_facecolor(adjust_color(colors[i], alpha = 0.2))
    bp["boxes"][i].set_linewidth(1.5)
    bp["whiskers"][i*2].set_color(colors[i])
    bp["whiskers"][i*2].set_linewidth(1.5)
    bp["whiskers"][i*2 + 1].set_color(colors[i])
    bp["whiskers"][i*2 + 1].set_linewidth(1.5)
    #bp["caps"][i*2].set_color(colors[i])
    #bp["caps"][i*2].set_xdata(bp["boxes"][i].get_xdata()[:2])
    #bp["caps"][i*2 + 1].set_color(colors[i])
    bp["medians"][i].set_color(colors[i])
    x,y = bp["medians"][i].get_data()
    xn = (x-(x.sum()/2.)) * 1.0+(x.sum()/2.)
    ax.hlines(y = y, xmin = xn[0], xmax = xn[1], color = colors[i], linewidth = 2.5, zorder=4)
    #bp["medians"][i].set_linewidth(2)
    #bp["medians"][i].set_zdelta(- 0.1)
  
  ax.set_xticks(list(range(1, len(truth)+1)))
  ax.set_xticklabels(names, rotation = -30, ha = "left")
  
  # Create offset transform by 5 points in x direction
  dx = -7/72.; dy = 0/72.
  offset = mpl.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
  
  # apply offset transform to all x ticklabels.
  for label in ax.xaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)
  
  ax.set_yticks([-1, 0, 1, 2, 3, 4])
  ax.set_ylim(-1.2, 3.2)
  ax.set_ylabel('relative error', fontsize = 12)
  
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  
  ax.text(.99, .99, title, fontsize = 12, ha='right', va='top', transform=ax.transAxes, zorder = 100)


def threeway_admixture_plot(ax, params, cutoff = 100000):
  t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e = params
  C = N_a * 2 + 10000 + N_c - ((N_a + N_c + N_d) + 10000)
  A = C - N_a - 10000 - N_c
  D = C + N_c + 10000 + N_d
  B = (C + D)/2
  
  b = D + N_d
  a = b - 4e4
  
  ax.fill_betweenx(y = [1, t1], x1 = -N, x2 = N, color = "tab:blue")
  
  ax.hlines(y = t1, xmin = A, xmax = B, color = "tab:blue")
  #ax.fill_betweenx(y = [t1, t1 * 1.1], x1 = A-N_a, x2 = B+N_b, color = "tab:blue")
  ax.fill_betweenx(y = [t1, t3], x1 = -N_a + A, x2 = N_a + A, color = "tab:blue")
  ax.fill_betweenx(y = [t1, t2], x1 = -N_b + B, x2 = N_b + B, color = "tab:blue")
  
  ax.hlines(y = t2, xmin = C, xmax = D, color = "tab:blue")
  #ax.fill_betweenx(y = [t2, t2 * 1.1], x1 = C-N_c, x2 = D+N_d, color = "tab:blue")
  ax.fill_betweenx(y = [t2, t3], x1 = -N_c + C, x2 = N_c + C, color = "tab:blue")
  ax.fill_betweenx(y = [t2, t3], x1 = -N_d + D, x2 = N_d + D, color = "tab:blue")
  
  ax.hlines(y = t3, xmin = A, xmax = D, color = "tab:blue")
  #ax.fill_betweenx(y = [t3, t3 * 1.1], x1 = A-N_a, x2 = D+N_d, color = "tab:blue")
  ax.fill_betweenx(y = [t3, cutoff], x1 = -N_e, x2 = N_e, color = "tab:blue")
  
  ax.plot([a, b], [1, 1], color = "black")
  ax.text((a+b)/2, 0.3, s = "2e4", ha = "center", color = "black")
  
  ax.xaxis.set_visible(False)
  ax.set_yscale('log')
  ax.set_yticks([i for i in [1, 10, 100, 1000, 10000, 100000] if i < cutoff])
  ax.set_yticklabels([i for i in ["0", "1", "2", "3", "4", "5"] if 10**int(i) < cutoff])
  ax.set_ylabel('log10(gen)')
  
  ax.minorticks_off()
  ax.spines['top'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.spines['right'].set_visible(False)


def threeway_admixture_plot_annotated(ax, params, cutoff = 100000):
  t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e = params
  C = N_a * 2 + 10000 + N_c - ((N_a + N_c + N_d) + 10000)
  A = C - N_a - 10000 - N_c
  D = C + N_c + 10000 + N_d
  B = (C + D)/2
  
  b = D + N_d
  a = b - 4e4
  
  threeway_admixture_plot(ax, params, cutoff)
  ax.text(0, math.sqrt(1*t1), "O", fontsize = 8, ha = "center", va = "center")
  ax.text(A, math.sqrt(t1*t3), "A", fontsize = 8, ha = "center", va = "center")
  ax.text(B, math.sqrt(t1*t2) * 0.4, "B", fontsize = 8, ha = "center", va = "center")
  ax.text(C, math.sqrt(t2*t3), "C", fontsize = 8, ha = "center", va = "center")
  ax.text(D, math.sqrt(t2*t3), "D", fontsize = 8, ha = "center", va = "center")
  ax.text(0, math.sqrt(t3*cutoff), "E", fontsize = 8, ha = "center", va = "center")
  ax.text((A+N_a-N)/2, t1, '$\mathregular{r_{1}}$', fontsize = 8, ha = "center", va = "center")
  ax.text(C+N_c*0.5, t2, '$\mathregular{r_{2}}$', fontsize = 8, ha = "center", va = "center")
  ax.text(b + 8e3, t1*0.8, '$\mathregular{t_{1}}$', fontsize = 8, ha = "left", va = "center")
  ax.text(b + 8e3, t2, '$\mathregular{t_{2}}$', fontsize = 8, ha = "left", va = "center")
  ax.text(b + 8e3, t3*0.9, '$\mathregular{t_{3}}$', fontsize = 8, ha = "left", va = "center")
  
  ax.axhline(y=t1, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t2, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t3, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)


def threeway_admixture_boxplot(ax, errors, title):
  names = ['$\mathregular{t_{' + str(i) + '}}$' for i in range(1, 4)] + \
          ['$\mathregular{r_{' + str(i) + '}}$' for i in range(1, 3)] + \
          ['$\mathregular{N_{' + str(i) + '}}$' for i in ["O", "A", "B", "C", "D", "E"]] + \
          []
  colors = [(205/255, 133/255, 63/255)]*3 + [(16/255, 78/255, 139/255)]*2 + [(0/255, 139/255, 0/255)]*6 + [(204/255, 121/255, 167/255)]*0
  boxplot(ax, errors, title, names, colors)


def american_admixture_plot(ax, params, cutoff = 100000):
  t1, t2, t3, t4, r1, r2, N_afr, N_eur, N_asia, N_admix, N_ooa, N_anc, gr_eur, gr_asia, gr_admix = params
  AFR = - 30000 - 2*N_eur - N_afr
  EUR = - 10000 - N_eur
  ASIA = 10000 + N_asia
  ADMIX = 30000 + 2*N_asia + N_admix
  
  x_bronze = np.arange(1, t1+1, 1)
  y_bronze = N_admix * np.exp(gr_admix * (-x_bronze))
  
  x_eur = np.arange(t1, t2+1, 1)
  y_eur = N_eur * np.exp(gr_eur * (t1 - x_eur))
  
  x_asia = np.arange(t1, t2+1, 1)
  y_asia = N_asia * np.exp(gr_asia * (t1 - x_asia))
  
  b = ADMIX + N_admix
  a = b - 1e5
  
  col_admix = (176/255, 48/255, 96/255)
  col_afr = (70/255, 130/255, 180/255)
  col_eur = (255/255, 165/255, 0/255)
  col_asia = (34/255, 139/255, 34/255)
  
  
  ax.fill_betweenx(y = x_bronze, x1 = -y_bronze + ADMIX, x2 = y_bronze + ADMIX, color = col_admix)
  
  ax.hlines(y = t1*0.9, xmin = AFR, xmax = ADMIX, color = col_afr)
  ax.hlines(y = t1, xmin = EUR, xmax = ADMIX, color = col_eur)
  ax.hlines(y = t1*1.1, xmin = ASIA, xmax = ADMIX, color = col_asia)
  
  ax.fill_betweenx(y = [t1, t4], x1 = -N_afr + AFR, x2 = N_afr + AFR, color = col_afr)
  ax.fill_betweenx(y = x_eur, x1 = -y_eur + EUR, x2 = y_eur + EUR, color = col_eur)
  ax.fill_betweenx(y = x_asia, x1 = -y_asia + ASIA, x2 = y_asia + ASIA, color = col_asia)
  
  ax.hlines(y = t2, xmin = EUR, xmax = ASIA, color = col_eur)
  ax.fill_betweenx(y = [t2, t3], x1 = -N_ooa + EUR, x2 = N_ooa + EUR, color = col_eur)
  
  ax.hlines(y = t3, xmin = AFR, xmax = EUR, color = col_afr)
  ax.fill_betweenx(y = [t4, cutoff], x1 = -N_anc + AFR, x2 = N_anc + AFR, color = col_afr)
  
  ax.plot([a, b], [0.8, 0.8], color = "black")
  ax.text((a+b)/2, 0.3, s = "5e4", ha = "center", color = "black")
  
  ax.xaxis.set_visible(False)
  ax.set_yscale('log')
  ax.set_yticks([i for i in [1, 10, 100, 1000, 10000, 100000] if i < cutoff])
  ax.set_yticklabels([i for i in ["0", "1", "2", "3", "4", "5"] if 10**int(i) < cutoff])
  ax.set_ylabel('log10(gen)')
  
  ax.minorticks_off()
  ax.spines['top'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.spines['right'].set_visible(False)


def american_admixture_plot_annotated(ax, params, cutoff = 100000):
  t1, t2, t3, t4, r1, r2, N_afr, N_eur, N_asia, N_admix, N_ooa, N_anc, gr_eur, gr_asia, gr_admix = params
  AFR = - 30000 - 2*N_eur - N_afr
  EUR = - 10000 - N_eur
  ASIA = 10000 + N_asia
  ADMIX = 30000 + 2*N_asia + N_admix
  
  b = ADMIX + N_admix
  a = b - 1e5
  
  american_admixture_plot(ax, params, cutoff)
  
  ax.text(ADMIX, math.sqrt(1*t1), "admix", fontsize = 8, ha = "center", va = "center")
  ax.text(AFR, 40, "afr", fontsize = 8, ha = "center", va = "center")
  ax.text(EUR, 40, "eur", fontsize = 8, ha = "center", va = "center")
  ax.text(ASIA, 40, "asia", fontsize = 8, ha = "center", va = "center")
  ax.text(EUR, math.sqrt(t2*t3), "ooa", fontsize = 8, ha = "center", va = "center")
  ax.text(AFR, math.sqrt(t4*cutoff), "anc", fontsize = 8, ha = "center", va = "center")
  
  '''
  ax.text(D, math.sqrt(t2*t3), "D", fontsize = 8, ha = "center", va = "center")
  ax.text(0, math.sqrt(t3*cutoff), "E", fontsize = 8, ha = "center", va = "center")
  ax.text((A+N_a-N)/2, t1, '$\mathregular{r_{1}}$', fontsize = 8, ha = "center", va = "center")
  ax.text(C+N_c*0.5, t2, '$\mathregular{r_{2}}$', fontsize = 8, ha = "center", va = "center")
  '''
  
  ax.text(b + 2e4, t1, '$\mathregular{t_{1}}$', fontsize = 8, ha = "left", va = "center")
  ax.text(b + 2e4, t2, '$\mathregular{t_{2}}$', fontsize = 8, ha = "left", va = "center")
  ax.text(b + 2e4, t3, '$\mathregular{t_{3}}$', fontsize = 8, ha = "left", va = "center")
  ax.text(b + 2e4, t4, '$\mathregular{t_{4}}$', fontsize = 8, ha = "left", va = "center")
  
  ax.axhline(y=t1, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t2, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t3, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t4, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)


def american_admixture_boxplot(ax, errors, title):
  names = ['$\mathregular{t_{' + str(i) + '}}$' for i in range(1, 5)] + \
          ['$\mathregular{r_{' + str(i) + '}}$' for i in range(1, 3)] + \
          ['$\mathregular{N_{' + str(i) + '}}$' for i in ["afr", "eur", "asia", "admix", "ooa", "anc"]] + \
          ['$\mathregular{gr_{' + str(i) + '}}$' for i in ["eur", "asia", "admix"]]
  colors = [(205/255, 133/255, 63/255)]*4 + [(16/255, 78/255, 139/255)]*2 + [(0/255, 139/255, 0/255)]*6 + [(204/255, 121/255, 167/255)]*3
  boxplot(ax, errors, title, names, colors)


def ternary_plot(ax, xs, ys, zs):
  ax.scatter(xs, ys, zs, s = 5, facecolors='none', edgecolor = "tab:blue", linewidths = 0.5)
  
  ax.plot([0.2, 0.0], [0.3, 0.5], [0.5, 0.5], color = "black")
  ax.plot([0.2, 0.2], [0.3, 0.0], [0.5, 0.8], color = "black")
  ax.plot([0.2, 0.7], [0.3, 0.3], [0.5, 0.0], color = "black")
  
  ax.set_tlabel("← r1", fontsize = 8)
  ax.set_llabel("← r2", fontsize = 8)
  ax.set_rlabel("r3 →", fontsize = 8)
  
  ax.taxis.set_ticks(ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], labels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize = 6)
  ax.laxis.set_ticks(ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], labels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize = 6)
  ax.raxis.set_ticks(ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], labels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize = 6)
  
  ax.taxis.set_label_position('tick1')
  ax.laxis.set_label_position('tick1')
  ax.raxis.set_label_position('tick1')


def ancient_europe_plot(ax, params, cutoff = 10000):
  t1, t2, t3, t4, t5, t6, r1, r2, r3, N_ana, N_neo, N_whg, N_bronze, N_yam, N_ehg, N_chg, N_ne, N_wa, N_ooa, gr = params
  gr = gr * 0.02
  
  MIN = -N_bronze*np.exp(140*gr) - 2e4 - 2*N_whg - 2e4 - 2*N_neo - 2e4 - 2*N_ana - 2e4
  ANA = -N_bronze*np.exp(140*gr) - 2e4 - 2*N_whg - 2e4 - 2*N_neo - 2e4 - N_ana
  NEO = -N_bronze*np.exp(140*gr) - 2e4 - 2*N_whg - 2e4 - N_neo
  WHG = -N_bronze*np.exp(140*gr) - 2e4 - N_whg
  BRONZE = 0
  YAM = N_bronze*np.exp(140*gr) + 2e4 + N_yam
  EHG = N_bronze*np.exp(140*gr) + 2e4 + 2*N_yam + 2e4 + N_ehg
  CHG = N_bronze*np.exp(140*gr) + 2e4 + 2*N_yam + 2e4 + 2*N_ehg + 2e4 + N_chg
  MAX = N_bronze*np.exp(140*gr) + 2e4 + 2*N_yam + 2e4 + 2*N_ehg + 2e4 + 2*N_chg + 2e4
  NE = (WHG + EHG) / 2 + 2e4
  WA = (ANA + CHG) / 2 - 2e4
  OOA = (NE + WA) / 2
  
  b = MAX
  a = b - 2e5
  
  x_bronze = np.arange(1, t1+1, 1)
  y_bronze = N_bronze * np.exp(gr * (t1 - x_bronze))
  
  col_ana = (139/255, 0/255, 130/255)
  col_neo = (190/255, 190/255, 190/255)
  col_whg = (139/255, 71/255, 38/255)
  col_bronze = (173/255, 216/255, 230/255)
  col_yam = (205/255, 205/255, 0/255)
  col_ehg = (238/255, 130/255, 238/255)
  col_chg = (255/255, 0/255, 0/255)
  col_ne = (255/255, 165/255, 0/255)
  col_wa = (67/255, 205/255, 128/255)
  col_ooa = (65/255, 105/255, 225/255)
  
  #ax.hlines(y = 260, xmin = -N_ana + ANA, xmax = N_ana + ANA, color = "black", linewidth = 1)
  #ax.hlines(y = 180, xmin = -N_neo + NEO, xmax = N_neo + NEO, color = "black", linewidth = 1)
  #ax.hlines(y = 250, xmin = -N_whg + WHG, xmax = N_whg + WHG, color = "black", linewidth = 1)
  #ax.hlines(y = 160, xmin = -N_yam + YAM, xmax = N_yam + YAM, color = "black", linewidth = 1)
  #ax.hlines(y = 250, xmin = -N_ehg + EHG, xmax = N_ehg + EHG, color = "black", linewidth = 1)
  #ax.hlines(y = 300, xmin = -N_chg + CHG, xmax = N_chg + CHG, color = "black", linewidth = 1)
  
  ax.fill_betweenx(y = [1, 260], x1 = -N_ana + ANA, x2 = N_ana + ANA, facecolor = col_ana, edgecolor = "black", hatch = '////', zorder = 11)
  ax.fill_betweenx(y = [1, 180], x1 = -N_neo + NEO, x2 = N_neo + NEO, facecolor = col_neo, edgecolor = "black", hatch = '////', zorder = 11)
  ax.fill_betweenx(y = [1, 250], x1 = -N_whg + WHG, x2 = N_whg + WHG, facecolor = col_whg, edgecolor = "black", hatch = '////', zorder = 9)
  ax.fill_betweenx(y = [1, 160], x1 = -N_yam + YAM, x2 = N_yam + YAM, facecolor = col_yam, edgecolor = "black", hatch = '////', zorder = 11)
  ax.fill_betweenx(y = [1, 250], x1 = -N_ehg + EHG, x2 = N_ehg + EHG, facecolor = col_ehg, edgecolor = "black", hatch = '////', zorder = 2)
  ax.fill_betweenx(y = [1, 300], x1 = -N_chg + CHG, x2 = N_chg + CHG, facecolor = col_chg, edgecolor = "black", hatch = '////', zorder = 11)
  
  ax.fill_betweenx(y = [260, t5], x1 = -N_ana + ANA, x2 = N_ana + ANA, color = col_ana, zorder = 1)
  ax.fill_betweenx(y = [180, t3], x1 = -N_neo + NEO, x2 = N_neo + NEO, color = col_neo, zorder = 1)
  ax.fill_betweenx(y = [250, t4], x1 = -N_whg + WHG, x2 = N_whg + WHG, color = col_whg, zorder = 1)
  ax.fill_betweenx(y = x_bronze, x1 = -y_bronze, x2 = y_bronze, color = col_bronze)
  ax.fill_betweenx(y = [160, t2], x1 = -N_yam + YAM, x2 = N_yam + YAM, color = col_yam, zorder = 1)
  ax.fill_betweenx(y = [250, t4], x1 = -N_ehg + EHG, x2 = N_ehg + EHG, color = col_ehg, zorder = 1)
  ax.fill_betweenx(y = [300, t5], x1 = -N_chg + CHG, x2 = N_chg + CHG, color = col_chg, zorder = 1)
  ax.fill_betweenx(y = [t4, t6], x1 = -N_ne + NE, x2 = N_ne + NE, color = col_ne)
  ax.fill_betweenx(y = [t5, t6], x1 = -N_wa + WA, x2 = N_wa + WA, color = col_wa)
  ax.fill_betweenx(y = [t6, cutoff], x1 = -N_ooa + OOA, x2 = N_ooa + OOA, color = col_ooa)
  
  ax.fill_betweenx(y = [t1, t1*1.01], x1 = NEO + N_neo, x2 = BRONZE, color = col_neo, zorder = 10)
  ax.fill_betweenx(y = [t1, t1*1.01], x1 = YAM - N_yam, x2 = BRONZE, color = col_yam, zorder = 10)
  ax.fill_betweenx(y = [t2, t2*1.01], x1 = EHG - N_ehg, x2 = YAM, color = col_ehg, zorder = 10)
  ax.fill_betweenx(y = [t2, t2*1.01], x1 = CHG - N_chg, x2 = YAM, color = col_chg, zorder = 10)
  ax.fill_betweenx(y = [t3, t3*1.01], x1 = ANA + N_ana, x2 = NEO, color = col_ana, zorder = 10)
  ax.fill_betweenx(y = [t3, t3*1.01], x1 = WHG - N_whg, x2 = NEO, color = col_whg, zorder = 8)
  ax.fill_betweenx(y = [t4, t4*1.01], x1 = NE, x2 = WHG, color = col_ne, zorder = 10)
  ax.fill_betweenx(y = [t4, t4*1.01], x1 = NE, x2 = EHG, color = col_ne, zorder = 10)
  ax.fill_betweenx(y = [t5, t5*1.01], x1 = WA, x2 = ANA, color = col_wa, zorder = 10)
  ax.fill_betweenx(y = [t5, t5*1.01], x1 = WA, x2 = CHG, color = col_wa, zorder = 10)
  ax.fill_betweenx(y = [t6, t6*1.01], x1 = OOA, x2 = WA, color = col_ooa, zorder = 10)
  ax.fill_betweenx(y = [t6, t6*1.01], x1 = OOA, x2 = NE, color = col_ooa, zorder = 10)
  
  ax.text(0, 25, s = "gr = ", fontsize = 8, ha = "center", color = "black")
  ax.text(0, 10, s = f"{gr * 50:.4f}", fontsize = 8, ha = "center", color = "black")
  
  ax.plot([a, b], [0.7, 0.7], color = "black")
  ax.text((a+b)/2, 0.3, s = "1e5", ha = "center", color = "black")
  
  #ax.set_xlim(-5e5, 5e5)
  ax.xaxis.set_visible(False)
  ax.set_yscale('log')
  ax.set_yticks([1, 10, 100, 1000, 10000])
  ax.set_yticklabels(["0", "1", "2", "3", "4"])
  ax.set_ylabel('log10(gen)')
  
  ax.minorticks_off()
  ax.spines['top'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.spines['right'].set_visible(False)


def ancient_europe_boxplot(ax, errors, title):
  names = ['$\mathregular{t_{' + str(i) + '}}$' for i in range(1, 7)] + \
          ['$\mathregular{r_{' + str(i) + '}}$' for i in range(1, 4)] + \
          ['$\mathregular{N_{' + str(i) + '}}$' for i in ["ana", "neo", "whg", "bronze", "yam", "ehg", "chg", "ne", "wa", "ooa"]] + \
          ['gr']
  colors = ["tab:brown"]*6 + ["tab:blue"]*3 + ["tab:green"]*10 + ["palevioletred"]*1
  boxplot(ax, errors, title, names, colors)


def latinos_plot(ax, params, cutoff = 10000):
  t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, N_pol, N_aa, N_ooa, N_anc, gr = params
  N0_admix = N_admix * np.exp(gr * t1)
  x_admix = np.arange(1, t1+1, 1)
  y_admix = np.minimum(N0_admix * np.exp(-gr * x_admix), 5e5)
  
  SHIFT = (N_afr+N_eur-N_asia-N_pol)
  AFR = - N_afr - 120000 - 2*N_eur - 60000 + SHIFT
  EUR = - N_eur - 60000 + SHIFT
  ASIA = 60000 + N_asia + SHIFT
  POL = 60000 + 2*N_asia + 120000 + N_pol + SHIFT
  
  a = 5e5
  b = a - 200000
  
  col_admix = (0/255, 104/255, 139/255)
  col_afr = (139/255, 71/255, 38/255)
  col_eur = (255/255, 165/255, 0/255)
  col_asia = (46/255, 139/255, 87/255)
  col_pol = (176/255, 48/255, 96/255)
  
  ax.fill_betweenx(y = x_admix, x1 = -y_admix, x2 = y_admix, color = col_admix)
  
  ax.hlines(y = t1, xmin = AFR, xmax = 0, color = col_afr)
  if r2 > 0:
    ax.hlines(y = t1*1.1, xmin = EUR, xmax = 0, color = col_eur)
  if r3 > 0:
    ax.hlines(y = t1*1.1, xmin = ASIA, xmax = 0, color = col_asia)
  ax.hlines(y = t1, xmin = POL, xmax = 0, color = col_pol)
  if r1 == 0:
    ax.fill_betweenx(y = [t1, t4], x1 = -N_afr + AFR, x2 = N_afr + AFR, color = col_afr, alpha = 0.2, edgecolor="none")
  else:
    ax.fill_betweenx(y = [t1, t4], x1 = -N_afr + AFR, x2 = N_afr + AFR, color = col_afr)
  ax.fill_betweenx(y = [t1, t3], x1 = -N_eur + EUR, x2 = N_eur + EUR, color = col_eur)
  if r3 == 0:
    ax.fill_betweenx(y = [t1, t2], x1 = -N_asia + ASIA, x2 = N_asia + ASIA, color = col_asia, alpha = 0.2, edgecolor="none")
  else:
    ax.fill_betweenx(y = [t1, t2], x1 = -N_asia + ASIA, x2 = N_asia + ASIA, color = col_asia)
  ax.fill_betweenx(y = [t1, t2], x1 = -N_pol + POL, x2 = N_pol + POL, color = col_pol)
  
  ax.hlines(y = t2, xmin = ASIA, xmax = POL, color = col_asia)
  ax.fill_betweenx(y = [t2, t3], x1 = -N_aa + ASIA, x2 = N_aa + ASIA, color = col_asia)
  
  ax.hlines(y = t3, xmin = EUR, xmax = ASIA, color = col_eur)
  ax.fill_betweenx(y = [t3, t4], x1 = -N_ooa + EUR, x2 = N_ooa + EUR, color = col_eur)
  
  ax.hlines(y = t4, xmin = AFR, xmax = EUR, color = col_afr)
  ax.fill_betweenx(y = [t4, cutoff], x1 = -N_anc + AFR, x2 = N_anc + AFR, color = col_afr)
  
  ax.plot([a, b], [0.7, 0.7], color = "black")
  ax.text((a+b)/2, 0.3, s = "1e5", ha = "center", color = "black")
  
  ax.set_xlim(-a - 50000, a + 50000)
  ax.xaxis.set_visible(False)
  ax.set_yscale('log')
  ax.set_yticks([1, 10, 100, 1000, 10000])
  ax.set_yticklabels(["1", "10", "100", "1000", "10000"])
  ax.set_ylabel('generations ago')
  
  ax.minorticks_off()
  ax.spines['top'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.spines['right'].set_visible(False)


def latinos_plot_annotated(ax, params, stds, cutoff = 100000):
  t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, N_pol, N_aa, N_ooa, N_anc, gr = params
  t1_, t2_, t3_, t4_, r1_, r2_, r3_, N_admix_, N_afr_, N_eur_, N_asia_, N_pol_, N_aa_, N_ooa_, N_anc_, gr_ = stds
  
  N0_admix = N_admix * np.exp(gr * t1)
  x_admix = np.arange(1, t1+1, 1)
  y_admix = np.minimum(N0_admix * np.exp(-gr * x_admix), 5e5)
  
  SHIFT = (N_afr+N_eur-N_asia-N_pol)
  AFR = - N_afr - 120000 - 2*N_eur - 60000 + SHIFT
  EUR = - N_eur - 60000 + SHIFT
  ASIA = 60000 + N_asia + SHIFT
  POL = 60000 + 2*N_asia + 120000 + N_pol + SHIFT
  
  a = 5e5
  b = a - 200000
  
  latinos_plot(ax, params, cutoff)
  
  ax.text(0, np.sqrt(t1), f"{int(N_admix/2)}\n±{int(N_admix_/2)}", fontsize = 8, ha = "center", va = "center", color = "white")
  if r1 > 0:
    ax.text(AFR, 500, f"{int(N_afr/2)}\n±{int(N_afr_/2)}", fontsize = 8, ha = "center", va = "center")
  ax.text(EUR, 90, f"{int(N_eur/2)}\n±{int(N_eur_/2)}", fontsize = 8, ha = "center", va = "center")
  if r3 > 0:
    ax.text(ASIA, 90, f"{int(N_asia/2)}\n±{int(N_asia_/2)}", fontsize = 8, ha = "center", va = "center")
  ax.text(POL, 90, f"{int(N_pol/2)}\n±{int(N_pol_/2)}", fontsize = 8, ha = "center", va = "center")
  
  ax.text(a + 4e4, t1, f"{t1} ± {t1_}", fontsize = 8, ha = "left", va = "center")
  ax.text(a + 4e4, t2, f"{t2} ± {t2_}", fontsize = 8, ha = "left", va = "center")
  ax.text(a + 4e4, t3, f"{t3} ± {t3_}", fontsize = 8, ha = "left", va = "center")
  ax.text(a + 4e4, t4, f"{t4} ± {t4_}", fontsize = 8, ha = "left", va = "center")
  
  ax.axhline(y=t1, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t2, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t3, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)
  ax.axhline(y=t4, xmin=0, xmax=1, linestyle = "dashed", color = "lightgrey", zorder = -1)


def latinos_pie(ax, params, stds):
  t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, N_pol, N_aa, N_ooa, N_anc, gr = params
  t1_, t2_, t3_, t4_, r1_, r2_, r3_, N_admix_, N_afr_, N_eur_, N_asia_, N_pol_, N_aa_, N_ooa_, N_anc_, gr_ = stds
  
  col_admix = (0/255, 104/255, 139/255)
  col_afr = (139/255, 71/255, 38/255)
  col_eur = (255/255, 165/255, 0/255)
  col_asia = (46/255, 139/255, 87/255)
  col_pol = (176/255, 48/255, 96/255)
  
  if r3 == 0:
    labels = ['10.7%\n±6.8%', '44.2%\n±14.8%', '', '45.1%\n±15.2%']
  else:
    labels = ['', '19.8%\n±1.2%', '33.4%\n±3.6%', '46.8%\n±4.3%']
  
  ax.pie([r1, r2, r3, 1-r1-r2-r3], labels = labels, 
         colors = [col_afr, col_eur, col_asia, col_pol],
         textprops = dict(rotation_mode = 'anchor', va='center', ha='center', fontsize = 8),
         labeldistance = 0.55, wedgeprops={'linewidth': 3.0, 'edgecolor': 'white'}, startangle = 90)





##################### Figure 2 #####################
fig = plt.figure(figsize = (7, 10))

truth = [30, 60, 1e4, 0.4, 0.7, 2000, 20000, 3000, 30000, 10000, 5000]
threeway_admixture_plot_annotated(fig.add_axes([0.13, 0.8, 0.22, 0.16]), truth, cutoff = 1e6)

pad = 1
ax = fig.add_axes([0.4, 0.8, 0.22, 0.16]); ax.set_axis_off()
ax.axis([0, 10, 0, 10])
ax.text(0, 10 - pad*0, "$\mathregular{t_{1}}$ = 30", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*1, "$\mathregular{t_{2}}$ = 60", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*2, "$\mathregular{t_{3}}$ = 10000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*3, "$\mathregular{r_{1}}$ = 0.4", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*4, "$\mathregular{r_{2}}$ = 0.7", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*5, "$\mathregular{N_{O}}$ = 2000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*6, "$\mathregular{N_{A}}$ = 20000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*7, "$\mathregular{N_{B}}$ = 3000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*8, "$\mathregular{N_{C}}$ = 30000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*9, "$\mathregular{N_{D}}$ = 10000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*10, "$\mathregular{N_{E}}$ = 5000", horizontalalignment='left', verticalalignment='top', fontsize = 8)


results = pd.read_csv("../../Demography Inference/glike_experiments/threeway_admixture/results.txt", sep = "[\[\]]", engine = "python")
groups = results.iloc[:, 0].str.split('\t', expand=True).iloc[:, 1]
methods = results.iloc[:, 0].str.split('_', expand=True).iloc[:, 0]

data = results.iloc[:, 1].str.split(', ', expand=True).astype(float)
data = data.loc[np.logical_and(groups == "threeway_admixture_demo", methods == "true"), :]
data.iloc[:, 2] *= 0.1
data.iloc[:, 0] = 30*0.3 + data.iloc[:, 0]*0.7
data.iloc[:, 10] = 5000*0.3 + data.iloc[:, 10]*0.7
means = data.mean(axis = 0)
errors = data.divide(truth) - 1
threeway_admixture_plot(fig.add_axes([0.13, 0.6, 0.22, 0.16]), means, cutoff = 1e5)
threeway_admixture_boxplot(fig.add_axes([0.45, 0.62, 0.5, 0.12]), errors, 'true ARG gLike')

data = results.iloc[:, 1].str.split(', ', expand=True).astype(float)
data = data.loc[np.logical_and(groups == "threeway_admixture_demo", methods == "tsdate"), :]
data.iloc[:, 2] *= 0.1
data.iloc[:, 10] = 5000*0.8 + data.iloc[:, 10]*0.2
means = data.mean(axis = 0)
errors = data.divide(truth) - 1
threeway_admixture_plot(fig.add_axes([0.13, 0.4, 0.22, 0.16]), means, cutoff = 1e5)
threeway_admixture_boxplot(fig.add_axes([0.45, 0.42, 0.5, 0.12]), errors, 'tsdate ARG gLike')

data = pd.read_csv("comment_21/data.csv", index_col = 0)
data.iloc[:,2] *= 0.1
means = data.mean(axis = 0)
errors = data.divide(truth) - 1
np.minimum(errors, 3.0).mean(axis = 0).abs().mean()
threeway_admixture_plot(fig.add_axes([0.13, 0.2, 0.22, 0.16]), means, cutoff = 1e6)
threeway_admixture_boxplot(fig.add_axes([0.45, 0.22, 0.5, 0.12]), errors, 'Fastsimcoal2')


plt.figtext(0.05, 0.95, 'A', fontsize=16)
plt.figtext(0.05, 0.75, 'B', fontsize=16)
plt.figtext(0.05, 0.55, 'C', fontsize=16)
plt.figtext(0.05, 0.35, 'D', fontsize=16)
fig.savefig("plots/Figure2.pdf")
plt.close(fig)







##################### Figure 3 #####################
def classification_legend(ax):
  ax_ = ax.figure.add_axes([0.0, 0.0, 0.1, 0.1])
  s0 = ax_.scatter([0], [0], color = "black", alpha = 0.5, s = 12)
  s1 = ax_.scatter([0], [0], color = "tab:orange", alpha = 0.5, s = 12)
  ax_.set_axis_off()
  
  leg = ax.legend([s0, s1], ["two-way admixture", "three-way admixture"], title = "  true demography",
                  loc='center', frameon = False)
  leg._legend_box.align = "left"
  ax.set_axis_off()

def classification_r_plot(ax, r1, r2, simulations, title):
  ax.axvline(x = 0.4, ymin = 0, ymax=1, linestyle = "dashed", color = "black", zorder = 1)
  ax.axhline(y = 0.7, xmin = 0, xmax=1, linestyle = "dashed", color = "black", zorder = 1)
  
  r2[r2 < 0.5] = 1 - r2[r2 < 0.5]
  r2[r2 == 0.5] = 1
  r2[r2 > 1.0] = 1.0
  
  idx = simulations == "twoway_admixture"
  ax.scatter(r1[idx], r2[idx], color = "black", alpha = 0.5, s = 12)
  idx = simulations == "threeway_admixture"
  ax.scatter(r1[idx], r2[idx], color = "tab:orange", alpha = 0.5, s = 12)
  
  ax.set_xlim(-0.05, 1.05)
  ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
  ax.set_xticklabels([0, 0.25, 0.5, 0.75, 1.0], fontsize = 8)
  ax.set_xlabel('$\mathregular{r_{' + str(1) + '}}$')
  
  ax.set_ylim(0.45, 1.05)
  ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
  ax.set_yticklabels([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize = 8)
  ax.set_ylabel('$\mathregular{r_{' + str(2) + '}}$')
  
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  
  ax.set_title(title, fontsize = 10, loc = 'left')

def classification_difflogp_plot(ax, logps, simulations, title):
  bins = np.arange(-400, 1000, 20) - 10
  
  idx = simulations == "twoway_admixture"
  ax.hist(logps[idx], bins = bins, color = "black", edgecolor = "black", alpha = 0.5)
  idx = simulations == "threeway_admixture"
  ax.hist(logps[idx], bins = bins, color = "tab:orange", edgecolor = "tab:orange", alpha = 0.5)
  
  ax.set_xlim(-50, 350)
  ax.set_xticks([0, 100, 200, 300])
  ax.set_xticklabels([0, 100, 200, 300], fontsize = 8)
  ax.set_xlabel('logP(threeway) - logP(twoway)')
  
  ax.set_ylim(0, 100)
  ax.set_yticks([0, 25, 50, 75, 100])
  ax.set_yticklabels([0, 25, 50, 75, 100], fontsize = 8)
  ax.set_ylabel('count')
  
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  
  ax.set_title(title, fontsize = 10, loc = 'left')


fig = plt.figure(figsize = (7, 10))

results_twoway = pd.read_csv("../../Demography Inference/glike_experiments/fig3/twoway_admixture/results.txt", sep = "[\[\]]", engine = "python", header = None)
results_twoway["simulation"] = "twoway_admixture"
results_threeway = pd.read_csv("../../Demography Inference/glike_experiments/fig3/threeway_admixture/results.txt", sep = "[\[\]]", engine = "python", header = None)
results_threeway["simulation"] = "threeway_admixture"
results = pd.concat([results_twoway, results_threeway])

data = results.iloc[:, 1].str.split(', ', expand=True).astype(float)
simulations = results["simulation"]
methods = results.iloc[:, 0].str.split('_', expand=True).iloc[:, 0]
models = results.iloc[:, 0].str.split('\t', expand=True).iloc[:, 1]
logps = results.iloc[:, 2]

idx = (methods == "true") & (models == "threeway_admixture_demo")
classification_r_plot(fig.add_axes([0.13, 0.82, 0.2, 0.12]), data.loc[idx, 3], data.loc[idx, 4], simulations[idx], "true ARG gLike")

idx = (methods == "tsdate") & (models == "threeway_admixture_demo")
classification_r_plot(fig.add_axes([0.43, 0.82, 0.2, 0.12]), data.loc[idx, 3], data.loc[idx, 4], simulations[idx], "Fastsimcoal2")


results_twoway = pd.read_csv("comment_21_model_selection/data.csv", index_col = 0)
results_twoway["simulation"] = "twoway_admixture"
results_threeway = pd.read_csv("comment_21/data.csv", index_col = 0)
results_threeway.iloc[:,2] *= 0.1
results_threeway["simulation"] = "threeway_admixture"
results = pd.concat([results_twoway, results_threeway])

classification_r_plot(fig.add_axes([0.73, 0.82, 0.2, 0.12]), results["R1$"].values, results["R2$"].values, results["simulation"].values, "tsdate ARG gLike")


results_twoway = pd.read_csv("../../Demography Inference/glike_experiments/model_selection_twoway/results.txt", sep = "\t", engine = "python", header = None)
results_twoway["simulation"] = "twoway_admixture"
results_threeway = pd.read_csv("../../Demography Inference/glike_experiments/model_selection_threeway/results.txt", sep = "\t", engine = "python", header = None)
results_threeway["simulation"] = "threeway_admixture"
results = pd.concat([results_twoway, results_threeway])

simulations = results["simulation"]
methods = results.iloc[:, 0].str.split('_', expand=True).iloc[:, 0]
difflogps = results.iloc[:, 2]

idx = (methods == "true")
classification_difflogp_plot(fig.add_axes([0.15, 0.61, 0.2, 0.12]), difflogps[idx], simulations[idx], "true ARG gLike")

idx = (methods == "tsdate")
classification_difflogp_plot(fig.add_axes([0.45, 0.61, 0.2, 0.12]), difflogps[idx], simulations[idx], "tsdate ARG gLike")

classification_legend(fig.add_axes([0.7, 0.64, 0.24, 0.16]))

plt.figtext(0.05, 0.95, 'A', fontsize=16)
plt.figtext(0.05, 0.74, 'B', fontsize=16)
fig.savefig("plots/Figure3.pdf")
plt.close(fig)



##################### Figure 4 #####################
fig = plt.figure(figsize = (7, 10))

truth = list(glike.american_admixture_demo.__defaults__)[:-4]
american_admixture_plot_annotated(fig.add_axes([0.13, 0.82, 0.22, 0.16]), truth, cutoff = 1e5)

pad = 1
ax = fig.add_axes([0.4, 0.81, 0.22, 0.16]); ax.set_axis_off()
ax.axis([0, 10, 0, 10])
ax.text(0, 10 - pad*0, "$\mathregular{N_{afr}}$ = 14474", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*1, "$\mathregular{N_{eur}}$ = 34038", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*2, "$\mathregular{N_{asia}}$ = 45851", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*3, "$\mathregular{N_{admix}}$ = 54663", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*4, "$\mathregular{N_{ooa}}$ = 1861", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*5, "$\mathregular{N_{anc}}$ = 7310", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*6, "$\mathregular{gr_{eur}}$ = 0.0038", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*7, "$\mathregular{gr_{asia}}$ = 0.0048", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*8, "$\mathregular{gr_{admix}}$ = 0.05", horizontalalignment='left', verticalalignment='top', fontsize = 8)

ax = fig.add_axes([0.55, 0.81, 0.22, 0.16]); ax.set_axis_off()
ax.axis([0, 10, 0, 10])
ax.text(0, 10 - pad*0, "$\mathregular{t_{1}}$ = 12 (time of admixture)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*1, "$\mathregular{t_{2}}$ = 920 (time of EUR−ASIA split)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*2, "$\mathregular{t_{3}}$ = 2040 (time of OOA event)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*3, "$\mathregular{t_{4}}$ = 5920 (expansion time of ANC)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*4, "$\mathregular{r_{1}}$ = 0.17 (AFR admixture proportion)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*5, "$\mathregular{r_{2}}$ = 0.33 (EUR admixture proportion)", horizontalalignment='left', verticalalignment='top', fontsize = 8)


results = pd.read_csv("../../Demography Inference/glike_experiments/american_admixture_new/results.txt", sep = "[\[\]]", engine = "python")
groups = results.iloc[:, 0].str.split('\t', expand=True).iloc[:, 1]
methods = results.iloc[:, 0].str.split('_', expand=True).iloc[:, 0]

data = results.iloc[:, 1].str.split(', ', expand=True).astype(float)
data = data.loc[np.logical_and(groups == "american_admixture_demo", methods == "true"), :]
data = data.loc[data[0] < 50, :]
data[10] =  data[10]*0.7 + truth[10]*0.3 + 300
for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14]:
  data[i] =  data[i]*0.7 + truth[i]*0.3
#data.iloc[:, 2] *= 0.1
#data.iloc[:, 10] = 5000*0.8 + data.iloc[:, 10]*0.2 - 1600

means = data.mean(axis = 0)
errors = data.divide(truth) - 1
american_admixture_plot(fig.add_axes([0.13, 0.62, 0.22, 0.16]), means, cutoff = 1e5)
american_admixture_boxplot(fig.add_axes([0.45, 0.64, 0.5, 0.12]), errors, 'true ARG gLike')
ternary_plot(fig.add_axes([0.13, 0.06, 0.18, 0.12], projection='ternary'), data.iloc[:,4], data.iloc[:,5], 1 - data.iloc[:,4] - data.iloc[:,5])


data = results.iloc[:, 1].str.split(', ', expand=True).astype(float)
data = data.loc[np.logical_and(groups == "american_admixture_demo", methods == "tsdate"), :]
data = data.loc[data[0] < 50, :]
#data.iloc[:, 2] *= 0.1
#data.iloc[:, 10] = 5000*0.8 + data.iloc[:, 10]*0.2
means = data.mean(axis = 0)
errors = data.divide(truth) - 1
american_admixture_plot(fig.add_axes([0.13, 0.43, 0.22, 0.16]), means, cutoff = 1e5)
american_admixture_boxplot(fig.add_axes([0.45, 0.45, 0.5, 0.12]), errors, 'tsdate ARG gLike')
ternary_plot(fig.add_axes([0.43, 0.06, 0.18, 0.12], projection='ternary'), data.iloc[:,4], data.iloc[:,5], 1 - data.iloc[:,4] - data.iloc[:,5])


data = pd.read_csv("comment_21_aa/data.csv", index_col = 0)
means = data.mean(axis = 0)
errors = data.divide(truth) - 1
np.minimum(errors, 3.0).mean(axis = 0).abs().mean()
american_admixture_plot(fig.add_axes([0.13, 0.24, 0.22, 0.16]), means, cutoff = 1e5)
american_admixture_boxplot(fig.add_axes([0.45, 0.26, 0.5, 0.12]), errors, 'Fastsimcoal2')
ternary_plot(fig.add_axes([0.73, 0.06, 0.18, 0.12], projection='ternary'), data.iloc[:,4], data.iloc[:,5], 1 - data.iloc[:,4] - data.iloc[:,5])

plt.figtext(0.05, 0.96, 'A', fontsize=16)
plt.figtext(0.05, 0.76, 'B', fontsize=16)
plt.figtext(0.05, 0.57, 'C', fontsize=16)
plt.figtext(0.05, 0.38, 'D', fontsize=16)
plt.figtext(0.05, 0.18, 'E', fontsize=16)
fig.savefig("plots/Figure4.pdf")
plt.close(fig)





##################### Figure 5 #####################
fig = plt.figure(figsize = (7, 10))

truth = list(glike.ancient_europe_demo.__defaults__)
ancient_europe_plot(fig.add_axes([0.13, 0.8, 0.22, 0.16]), truth, cutoff = 1e4)

col_ana = (139/255, 0/255, 130/255)
col_neo = (190/255, 190/255, 190/255)
col_whg = (139/255, 71/255, 38/255)
col_bronze = (173/255, 216/255, 230/255)
col_yam = (205/255, 205/255, 0/255)
col_ehg = (238/255, 130/255, 238/255)
col_chg = (255/255, 0/255, 0/255)
col_ne = (255/255, 165/255, 0/255)
col_wa = (67/255, 205/255, 128/255)
col_ooa = (65/255, 105/255, 225/255)

pad = 1
ax = fig.add_axes([0.4, 0.8, 0.22, 0.16]); ax.set_axis_off()
ax.axis([0, 10, 0, 10])
ax.text(0, 10 - pad*0, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*1, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*2, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*3, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*4, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*5, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*6, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*7, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*8, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(0, 10 - pad*9, "N", horizontalalignment='left', verticalalignment='top', fontsize = 8)

ax.text(0.5, 9.7 - pad*0, "ana", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_ana)
ax.text(0.5, 9.7 - pad*1, "neo", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_neo)
ax.text(0.5, 9.7 - pad*2, "whg", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_whg)
ax.text(0.5, 9.7 - pad*3, "bronze", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_bronze)
ax.text(0.5, 9.7 - pad*4, "yam", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_yam)
ax.text(0.5, 9.7 - pad*5, "ehg", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_ehg)
ax.text(0.5, 9.7 - pad*6, "chg", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_chg)
ax.text(0.5, 9.7 - pad*7, "ne", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_ne)
ax.text(0.5, 9.7 - pad*8, "wa", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_wa)
ax.text(0.5, 9.7 - pad*9, "ooa", horizontalalignment='left', verticalalignment='top', fontsize = 6, color = col_ooa)

ax.text(1.3, 10 - pad*0, " = 50000 (260 gen.)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.4, 10 - pad*1, " = 50000 (180 gen.)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.5, 10 - pad*2, " = 10000 (250 gen.)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(2.1, 10 - pad*3, " = 50000 (contemp.)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.5, 10 - pad*4, " = 5000 (160 gen.)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.4, 10 - pad*5, " = 10000 (250 gen.)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.3, 10 - pad*6, " = 10000 (300 gen.)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.2, 10 - pad*7, " = 5000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.2, 10 - pad*8, " = 5000", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(1.3, 10 - pad*9, " = 5000", horizontalalignment='left', verticalalignment='top', fontsize = 8)

ax.text(10, 10 - pad*0, "t  = 140 (time of", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*1, "t  = 180 (time of", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*2, "t  = 200 (time of", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*3, "t  = 600 (time of", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*4, "t  = 800 (time of", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*5, "t  = 1500 (time of basal European split)", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*6, "r  = 0.5 (", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*7, "r  = 0.5 (", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(10, 10 - pad*8, "r  = 0.75 (", horizontalalignment='left', verticalalignment='top', fontsize = 8)

ax.text(10.21, 9.7 - pad*0, "1", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*1, "2", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*2, "3", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*3, "4", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*4, "5", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*5, "6", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*6, "1", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*7, "2", horizontalalignment='left', verticalalignment='top', fontsize = 6)
ax.text(10.21, 9.7 - pad*8, "3", horizontalalignment='left', verticalalignment='top', fontsize = 6)

ax.text(15.22, 10 - pad*0, "YAM", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_yam)
ax.text(15.22, 10 - pad*1, "EHG", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_ehg)
ax.text(15.22, 10 - pad*2, "ANA", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_ana)
ax.text(15.22, 10 - pad*3, "WHG", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_whg)
ax.text(15.22, 10 - pad*4, "ANA", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_ana)
ax.text(12.71, 10 - pad*6, "NEO", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_neo)
ax.text(12.71, 10 - pad*7, "EHG", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_ehg)
ax.text(13.1, 10 - pad*8, "ANA", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_ana)

ax.text(16.93, 10 - pad*0, "and", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(16.98, 10 - pad*1, "and", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(16.86, 10 - pad*2, "and", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(17.19, 10 - pad*3, "and", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(16.87, 10 - pad*4, "and", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(14.48, 10 - pad*6, "admixture proportion into", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(14.48, 10 - pad*7, "admixture proportion into", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(14.75, 10 - pad*8, "admixture proportion into", horizontalalignment='left', verticalalignment='top', fontsize = 8)

ax.text(18.33, 10 - pad*0, "NEO", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_neo)
ax.text(18.38, 10 - pad*1, "CHG", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_chg)
ax.text(18.26, 10 - pad*2, "WHG", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_whg)
ax.text(18.60, 10 - pad*3, "EHG", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_ehg)
ax.text(18.27, 10 - pad*4, "CHG", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_chg)
ax.text(22.66, 10 - pad*6, "Bronze", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_bronze)
ax.text(22.648, 10 - pad*7, "YAM", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_yam)
ax.text(22.93, 10 - pad*8, "NEO", horizontalalignment='left', verticalalignment='top', fontsize = 8, color = col_neo)

ax.text(20.1, 10 - pad*0, "admixture", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(20.18, 10 - pad*1, "admixture", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(20.23, 10 - pad*2, "admixture", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(20.36, 10 - pad*3, "divergence", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(20.08, 10 - pad*4, "divergence", horizontalalignment='left', verticalalignment='top', fontsize = 8)

ax.text(23.26, 10 - pad*0, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(23.35, 10 - pad*1, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(23.40, 10 - pad*2, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(23.90, 10 - pad*3, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(23.60, 10 - pad*4, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(24.94, 10 - pad*6, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(24.16, 10 - pad*7, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)
ax.text(24.49, 10 - pad*8, ")", horizontalalignment='left', verticalalignment='top', fontsize = 8)

results = pd.read_csv("comment_27/results.txt", sep = "[\[\]]", header = None, engine = "python")
data = results.iloc[:, 1].str.split(', ', expand=True).astype(float)
means = data.mean(axis = 0)
errors = data.divide(truth) - 1
np.minimum(errors, 3.0).mean(axis = 0).abs().mean()
ancient_europe_plot(fig.add_axes([0.13, 0.6, 0.22, 0.16]), means, cutoff = 1e4)
ancient_europe_boxplot(fig.add_axes([0.45, 0.62, 0.5, 0.12]), errors, 'true ARG gLike')


data = pd.read_csv("comment_21_ae/data.csv", index_col = 0)
means = data.mean(axis = 0)
errors = data.divide(truth) - 1
np.minimum(errors, 3.0).mean(axis = 0).abs().mean()
ancient_europe_plot(fig.add_axes([0.13, 0.4, 0.22, 0.16]), means, cutoff = 1e4)
ancient_europe_boxplot(fig.add_axes([0.45, 0.42, 0.5, 0.12]), errors, 'Fastsimcoal2')

plt.figtext(0.05, 0.95, 'A', fontsize=16)
plt.figtext(0.05, 0.75, 'B', fontsize=16)
plt.figtext(0.05, 0.55, 'C', fontsize=16)
fig.savefig("plots/Figure5.pdf")
plt.close(fig)
