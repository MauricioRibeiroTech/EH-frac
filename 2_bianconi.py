import matplotlib.pyplot as plt
import pandas as pd

# Dados extraídos dos seus arquivos
data = {
    'N': [100, 200, 500, 1000],
    'Mech': [8.264473728817745, 9.702674918258861, 11.576196169404573, 12.979399059721818],
    'Elec': [7.179534101548926, 9.063693796020573, 11.344878667460273, 12.535029426546338]
}
df = pd.DataFrame(data)

# Configurações Estéticas
plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'axes.linewidth': 0.8})
fig, ax = plt.subplots(figsize=(4.5, 3.5))

# Plotagem
ax.plot(df['N'], df['Mech'], 'o-', label='Mechanical Layer (Small-World)', color='#2E4053', lw=1.5, markersize=6)
ax.plot(df['N'], df['Elec'], 's--', label='Electrical Layer (Random)', color='#C0392B', lw=1.5, markersize=6)

# Estilização
ax.set_xscale('log') # Escala logarítmica é padrão para N em redes
ax.set_xlabel('Network Size ($N$)')
ax.set_ylabel('Structural Entropy ($\Sigma$)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, which="both", ls="--", alpha=0.3)
ax.legend(frameon=False)

plt.tight_layout()
plt.savefig('bianconi_entropy_scaling.png', dpi=600)
plt.show()
