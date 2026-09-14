import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Cargar el archivo de resultados
# ---------------------------------------------------------
if len(sys.argv) < 2:
    sys.exit("Uso: python plot_benchmark.py <ruta_al_archivo_cell_benchmark_results.txt>")

results_path = sys.argv[1]

if not os.path.exists(results_path):
    sys.exit(f"Error: El archivo '{results_path}' no existe.")

# Cargar CSV generado por el script de rendimiento
df = pd.read_csv(results_path)

# Asegurar orden por número de puntos
df = df.sort_values(by="TotalPoints")

# ---------------------------------------------------------
# Crear directorio de salida para las imágenes
# ---------------------------------------------------------
output_dir = os.path.dirname(results_path)
plot_filename = os.path.join(output_dir, "benchmark_performance_plots.png")

# ---------------------------------------------------------
# Generar Gráficos
# ---------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
fig.suptitle("Análisis de Rendimiento - Voronoi Mesh Scaling", fontsize=16, fontweight='bold')

# Obtener mallas únicas si hay varias
meshes = df["MeshName"].unique()

for mesh in meshes:
    mesh_df = df[df["MeshName"] == mesh]
    
    # --- Gráfico 1: Puntos vs Tiempo Total ---
    axes[0].plot(
        mesh_df["TotalPoints"], 
        mesh_df["TotalTime"], 
        marker='o', 
        linewidth=2, 
        label=f"{mesh} (Total Time)"
    )
    
    # --- Gráfico 2: Desglose por Etapas (Puntos vs Voronoi) ---
    axes[1].plot(
        mesh_df["TotalPoints"], 
        mesh_df["TimePointGen"], 
        marker='s', 
        linestyle='--', 
        label=f"{mesh} (Point Gen)"
    )
    axes[1].plot(
        mesh_df["TotalPoints"], 
        mesh_df["TimeVoronoiGen"], 
        marker='^', 
        linestyle='-', 
        label=f"{mesh} (Voronoi Gen)"
    )

# Configuración Gráfico 1
axes[0].set_title("Tiempo Total vs Número de Puntos", fontsize=12)
axes[0].set_xlabel("Número de Puntos (TotalPoints)", fontsize=11)
axes[0].set_ylabel("Tiempo de Ejecución (Segundos)", fontsize=11)
axes[0].grid(True, linestyle='--', alpha=0.6)
axes[0].legend()

# Configuración Gráfico 2
axes[1].set_title("Desglose de Tiempos por Etapa", fontsize=12)
axes[1].set_xlabel("Número de Puntos (TotalPoints)", fontsize=11)
axes[1].grid(True, linestyle='--', alpha=0.6)
axes[1].legend()

plt.tight_layout()

# Guardar y mostrar
plt.savefig(plot_filename, dpi=300)
print(f"Gráfico guardado exitosamente en: {plot_filename}")
plt.show()