import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

if len(sys.argv) < 2:
    sys.exit("Usage: python plotBenchmark.py <path_to_benchmark_results.txt>")

file_path = sys.argv[1]

if not os.path.exists(file_path):
    sys.exit(f"Error: File '{file_path}' does not exist.")

# Load benchmark dataset
df = pd.read_csv(file_path)

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
fig, axes = plt.subplots(2, 2, figsize=(15, 11))
fig.suptitle('Dask Voronoi Performance & Scaling Benchmark', fontsize=16, fontweight='bold')

# Plot 1: Total Execution Time vs Processors
sns.lineplot(
    data=df, x='NProc', y='TotalTime', hue='MeshName', 
    marker='o', linewidth=2.5, ax=axes[0, 0]
)
axes[0, 0].set_title('Total Execution Time vs Processors', fontsize=12, fontweight='bold')
axes[0, 0].set_xlabel('Number of Processors (NProc)')
axes[0, 0].set_ylabel('Total Time (seconds)')
axes[0, 0].set_xticks(range(1, df['NProc'].max() + 1))

# Plot 2: Speedup Factor vs Ideal Scaling
sns.lineplot(
    data=df, x='NProc', y='Speedup', hue='MeshName', 
    marker='s', linewidth=2.5, ax=axes[0, 1]
)
max_nproc = df['NProc'].max()
axes[0, 1].plot([1, max_nproc], [1, max_nproc], 'k--', label='Ideal Speedup', alpha=0.7)
axes[0, 1].set_title('Speedup Factor vs Processors', fontsize=12, fontweight='bold')
axes[0, 1].set_xlabel('Number of Processors (NProc)')
axes[0, 1].set_ylabel('Speedup (x Times Faster)')
axes[0, 1].set_xticks(range(1, max_nproc + 1))
axes[0, 1].legend()

# Plot 3: Time Breakdown per Stage
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

sns.barplot(
    data=df_melted, x='NProc', y='Time', hue='Stage', 
    errorbar=None, ax=axes[1, 0]
)
axes[1, 0].set_title('Execution Time Breakdown per Stage', fontsize=12, fontweight='bold')
axes[1, 0].set_xlabel('Number of Processors (NProc)')
axes[1, 0].set_ylabel('Time (seconds)')

# Plot 4: Speedup Heatmap Matrix
pivot_sp = df.pivot(index='MeshName', columns='NProc', values='Speedup')

# Ensure heatmap data is clean and valid
if not pivot_sp.empty and not pivot_sp.isna().all().all():
    sns.heatmap(pivot_sp, annot=True, fmt=".2f", cmap="YlGnBu", ax=axes[1, 1], cbar_kws={'label': 'Speedup Factor'})
    axes[1, 1].set_title('Speedup Heatmap Matrix', fontsize=12, fontweight='bold')
    axes[1, 1].set_xlabel('Number of Processors (NProc)')
    axes[1, 1].set_ylabel('Mesh Case')
else:
    axes[1, 1].text(0.5, 0.5, 'Insufficient data for heatmap', ha='center', va='center')

plt.tight_layout()

# Save image
output_image = os.path.join(os.path.dirname(file_path), "benchmark_performance.png")
plt.savefig(output_image, dpi=300)
print(f"Chart successfully saved to: {output_image}")
plt.show()