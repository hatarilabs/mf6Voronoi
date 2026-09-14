import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

if len(sys.argv) < 2:
    sys.exit("Usage: python caseNproc_plotBenchmark.py <path_to_benchmark_results.txt>")

file_path = sys.argv[1]

if not os.path.exists(file_path):
    sys.exit(f"Error: File '{file_path}' does not exist.")

# Load benchmark dataset
df = pd.read_csv(file_path)
print(df.head())

# Ensure correct data types to prevent merge errors
df['NProc'] = df['NProc'].astype(int)
df['TotalTime'] = df['TotalTime'].astype(float)

# Compute Speedup relative to normalRun (nproc=1)
speedup_rows = []
for mesh in df['MeshName'].unique():
    m_df = df[df['MeshName'] == mesh]
    base_time_row = m_df[m_df['NProc'] == 1]
    
    if not base_time_row.empty:
        base_time = base_time_row['TotalTime'].values[0]
        for _, row in m_df.iterrows():
            speedup = base_time / row['TotalTime'] if row['TotalTime'] > 0 else 0.0
            speedup_rows.append((row['MeshName'], int(row['NProc']), speedup))

df_speedup = pd.DataFrame(speedup_rows, columns=['MeshName', 'NProc', 'Speedup'])

# Merge speedup with base df
df = pd.merge(df, df_speedup, on=['MeshName', 'NProc'])

# Verify that df is not empty after merge
if df.empty:
    sys.exit("Error: The merged DataFrame is empty. Check if NProc=1 exists in your benchmark output file.")

# Style configuration
sns.set_theme(style="whitegrid", palette="deep")
output_dir = os.path.dirname(file_path)

# ---------------------------------------------------------
# Plot 1: Total Execution Time vs Processors
# ---------------------------------------------------------
plt.figure(figsize=(8, 6))
sns.lineplot(
    data=df, x='NProc', y='TotalTime', hue='MeshName', 
    marker='o', linewidth=2.5
)
plt.title('Total Execution Time vs Processors', fontsize=12, fontweight='bold')
plt.xlabel('Number of Processors (NProc)')
plt.ylabel('Total Time (seconds)')
plt.xticks(range(1, df['NProc'].max() + 1))
plt.tight_layout()

out_p1 = os.path.join(output_dir, "plot1_total_execution_time.png")
plt.savefig(out_p1, dpi=300)
plt.close()
print(f"Plot 1 saved to: {out_p1}")

# ---------------------------------------------------------
# Plot 2: Speedup Factor vs Ideal Scaling
# ---------------------------------------------------------
plt.figure(figsize=(8, 6))
sns.lineplot(
    data=df, x='NProc', y='Speedup', hue='MeshName', 
    marker='s', linewidth=2.5
)
max_nproc = df['NProc'].max()
plt.plot([1, max_nproc], [1, max_nproc], 'k--', label='Ideal Speedup', alpha=0.7)
plt.title('Speedup Factor vs Processors', fontsize=12, fontweight='bold')
plt.xlabel('Number of Processors (NProc)')
plt.ylabel('Speedup (x Times Faster)')
plt.xticks(range(1, max_nproc + 1))
plt.legend()
plt.tight_layout()

out_p2 = os.path.join(output_dir, "plot2_speedup_factor.png")
plt.savefig(out_p2, dpi=300)
plt.close()
print(f"Plot 2 saved to: {out_p2}")

# ---------------------------------------------------------
# Plot 3: Time Breakdown per Stage
# ---------------------------------------------------------
df_melted = df.melt(
    id_vars=['MeshName', 'NProc'], 
    value_vars=['TimePointGen', 'TimeVoronoiGen', 'TimeShapefileGen'],
    var_name='Stage', value_name='Time'
)
df_melted['Stage'] = df_melted['Stage'].map({
    'TimePointGen': 'Point Generation',
    'TimeVoronoiGen': 'Voronoi Generation',
    'TimeShapefileGen': 'Shapefile Generation'
})

plt.figure(figsize=(8, 6))
sns.barplot(
    data=df_melted, x='NProc', y='Time', hue='Stage', 
    errorbar=None
)
plt.title('Execution Time Breakdown per Stage', fontsize=12, fontweight='bold')
plt.xlabel('Number of Processors (NProc)')
plt.ylabel('Time (seconds)')
plt.tight_layout()

out_p3 = os.path.join(output_dir, "plot3_time_breakdown.png")
plt.savefig(out_p3, dpi=300)
plt.close()
print(f"Plot 3 saved to: {out_p3}")

# ---------------------------------------------------------
# Plot 4: Speedup Heatmap Matrix
# ---------------------------------------------------------
pivot_sp = df.pivot(index='MeshName', columns='NProc', values='Speedup')

plt.figure(figsize=(8, 6))
if not pivot_sp.empty and not pivot_sp.isna().all().all():
    sns.heatmap(pivot_sp, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Speedup Factor'})
    plt.title('Speedup Heatmap Matrix', fontsize=12, fontweight='bold')
    plt.xlabel('Number of Processors (NProc)')
    plt.ylabel('Mesh Case')
else:
    plt.text(0.5, 0.5, 'Insufficient data for heatmap', ha='center', va='center')

plt.tight_layout()

out_p4 = os.path.join(output_dir, "plot4_speedup_heatmap.png")
plt.savefig(out_p4, dpi=300)
plt.close()
print(f"Plot 4 saved to: {out_p4}")