import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams

# ==========================================
# 1. CONFIGURAÇÕES GLOBAIS - PADRÃO NATURE
# ==========================================
plt.style.use('default')
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 12
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 9
rcParams['axes.linewidth'] = 0.9 
rcParams['savefig.dpi'] = 600

def apply_clean_layout(ax):
    """Remove bordas superior/direita e adiciona grid sutil."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, linestyle='--', alpha=0.3, zorder=0)

# ==========================================
# 2. PARÂMETROS, CORES E SHAPES
# ==========================================
alphas = [0.6, 0.8, 1.0]
N_values = [100, 200, 500, 1000]

color_map = {100: '#1f77b4', 200: '#ff7f0e', 500: '#2ca02c', 1000: '#d62728'}
marker_map = {100: 'o', 200: 'v', 500: 's', 1000: 'D'}

metrics = [
    ('Kuramoto_R', r'Synchrony ($R$)'),
    ('Power_Avg', r'Harvested Power ($P_{avg}$)'),
    ('Shannon_H', r'Shannon Entropy ($H$)'),
    ('Mutual_Info', r'Mutual Information ($MI$)')
]
labels_abcd = ['(a)', '(b)', '(c)', '(d)']

# ==========================================
# 3. GERAÇÃO DOS PAINÉIS 2x2 (PARA CADA ALPHA)
# ==========================================
for alpha in alphas:
    datasets = {}
    for N in N_values:
        file_path = f"dados_alpha_{alpha}_N_{N}.txt"
        if os.path.exists(file_path):
            datasets[N] = pd.read_csv(file_path, sep='\t')
    
    if not datasets:
        continue

    fig, axs = plt.subplots(2, 2, figsize=(7.5, 6.0), sharex=True)
    axs = axs.flatten()

    for idx, (metric_key, metric_label) in enumerate(metrics):
        ax = axs[idx]
        apply_clean_layout(ax)

        for N in N_values:
            if N in datasets:
                df = datasets[N]
                y_data = df[metric_key].abs() if 'Shannon' in metric_key else df[metric_key]
                x_data = df['Sigma_Inter']
                
                step = max(1, len(df) // 15) 
                
                # Plot da linha contínua
                ax.plot(x_data, y_data, color=color_map[N], lw=1.5, alpha=0.8, zorder=2)
                
                # Plot dos shapes por cima
                ax.plot(x_data[::step], y_data[::step], label=f'N = {N}', 
                        linestyle='None', marker=marker_map[N], markersize=6, 
                        markerfacecolor=color_map[N], markeredgecolor='white', 
                        markeredgewidth=0.8, zorder=3)

        ax.set_ylabel(metric_label)
        
        # ---> LIMITES DINÂMICOS PARA R <---
        if metric_key == 'Kuramoto_R':
            ymin, ymax = ax.get_ylim()
            # Evita o offset 1e-8 nos painéis de alta sincronia, 
            # mas permite mostrar a curva descendo caso o R caia de verdade.
            ax.set_ylim(min(ymin, 0.95), max(ymax, 1.02))
        
        # Labels (a), (b), (c), (d)
        ax.text(-0.16, 1.05, labels_abcd[idx], transform=ax.transAxes, 
                fontsize=13, fontweight='bold', va='bottom', ha='right')

        if idx == 0:
            ax.legend(frameon=False, loc='best', handlelength=1.0)
        if idx >= 2:
            ax.set_xlabel(r'Piezoelectric Coupling ($\sigma_{inter}$)')

    plt.tight_layout()
    filename_panel = f"painel_convergencia_alpha_{alpha}.png"
    plt.savefig(filename_panel, bbox_inches='tight')
    plt.close(fig)
    print(f"Gerado: {filename_panel}")


# ==========================================
# 4. GRÁFICO DE DESTAQUE: COMPARAÇÃO GLOBAL DOS ALPHAS (N=1000)
# ==========================================
fig_iso, ax_iso = plt.subplots(figsize=(5.5, 4.0))
apply_clean_layout(ax_iso)

alpha_colors = {0.6: '#4B0082', 0.8: '#FF8C00', 1.0: '#008080'}
alpha_markers = {0.6: 'o', 0.8: 's', 1.0: '^'}
N_target = 1000

for alpha in alphas:
    file_path = f"dados_alpha_{alpha}_N_{N_target}.txt"
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t')
        x_data = df['Sigma_Inter']
        y_data = df['Power_Avg']
        step = max(1, len(df) // 12)
        
        ax_iso.plot(x_data, y_data, color=alpha_colors[alpha], lw=2.0, alpha=0.8, zorder=2)
        ax_iso.plot(x_data[::step], y_data[::step], label=rf'$\alpha = {alpha}$', 
                    linestyle='None', marker=alpha_markers[alpha], markersize=7, 
                    markerfacecolor=alpha_colors[alpha], markeredgecolor='white', 
                    markeredgewidth=1.0, zorder=3)

ax_iso.set_ylabel(r'Harvested Power ($P_{avg}$)', fontweight='bold')
ax_iso.set_xlabel(r'Piezoelectric Coupling ($\sigma_{inter}$)', fontweight='bold')
ax_iso.legend(frameon=False, loc='best', title=f'N = {N_target}', title_fontsize=11)

ax_iso.text(-0.16, 1.05, '(e)', transform=ax_iso.transAxes, 
            fontsize=13, fontweight='bold', va='bottom', ha='right')

plt.tight_layout()
filename_iso = "destaque_impacto_alpha_potencia.png"
plt.savefig(filename_iso, bbox_inches='tight')
plt.close(fig_iso)
print(f"Gerado: {filename_iso}")