#!/usr/bin/env python3
"""
Haplogroup Frequency Table Generator
Processes haplogrep3 output CSV and creates a frequency table similar to mitomap design

Output formats:
- txt: Plain text table with three-column layout
- html: Professional HTML table with CSS styling, color-coding, and responsive design
- both: Generates both text and HTML versions

Usage:
    python haplogroup_frequency_table.py input.csv -f html -o output_name
    python haplogroup_frequency_table.py input.csv -f both
    python haplogroup_frequency_table.py input.csv  # defaults to text output
"""

import pandas as pd
import re
from collections import defaultdict
import numpy as np

def extract_major_haplogroup(haplogroup_str):
    """
    Extract the major haplogroup from detailed haplogroup string
    Examples: H1a1a -> H, L3e1a -> L3, M7b1a -> M, etc.
    """
    if pd.isna(haplogroup_str) or haplogroup_str == '':
        return 'Unknown'
    
    # Remove any spaces and convert to uppercase
    hg = str(haplogroup_str).strip().upper()
    
    # Handle special cases and extract major haplogroup
    if hg.startswith('L'):
        # For L haplogroups, extract L + number (L0, L1, L2, etc.)
        match = re.match(r'(L\d+)', hg)
        if match:
            return match.group(1)
        else:
            return 'L'
    elif hg.startswith('HV'):
        return 'HV'
    elif hg.startswith('H'):
        return 'H'
    elif hg.startswith('U'):
        return 'U'
    elif hg.startswith('J'):
        return 'J'
    elif hg.startswith('T'):
        return 'T'
    elif hg.startswith('K'):
        return 'K'
    elif hg.startswith('V'):
        return 'V'
    elif hg.startswith('X'):
        return 'X'
    elif hg.startswith('W'):
        return 'W'
    elif hg.startswith('I'):
        return 'I'
    elif hg.startswith('N'):
        return 'N'
    elif hg.startswith('A'):
        return 'A'
    elif hg.startswith('B'):
        return 'B'
    elif hg.startswith('F'):
        return 'F'
    elif hg.startswith('R'):
        return 'R'
    elif hg.startswith('P'):
        return 'P'
    elif hg.startswith('Y'):
        return 'Y'
    elif hg.startswith('S'):
        return 'S'
    elif hg.startswith('O'):
        return 'O'
    elif hg.startswith('M'):
        return 'M'
    elif hg.startswith('D'):
        return 'D'
    elif hg.startswith('C'):
        return 'C'
    elif hg.startswith('G'):
        return 'G'
    elif hg.startswith('E'):
        return 'E'
    elif hg.startswith('Q'):
        return 'Q'
    elif hg.startswith('Z'):
        return 'Z'
    else:
        # For any other cases, try to extract the first letter(s)
        match = re.match(r'([A-Z]+)', hg)
        if match:
            return match.group(1)
        else:
            return 'Unknown'

def categorize_haplogroup(major_hg):
    """
    Categorize haplogroups into L (African), M (Asian), or N (Eurasian) lineages
    """
    # L lineages (African)
    l_lineages = {'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9'}
    
    # M lineages (Asian)
    m_lineages = {'M', 'D', 'C', 'G', 'E', 'Q', 'Z'}
    
    # N lineages (Eurasian) - everything else that's not L or M
    n_lineages = {'H', 'U', 'B', 'J', 'T', 'K', 'F', 'A', 'R', 'HV', 'N', 'I', 'V', 'X', 'W', 'P', 'Y', 'S', 'O'}
    
    if major_hg in l_lineages:
        return 'L'
    elif major_hg in m_lineages:
        return 'M'
    elif major_hg in n_lineages:
        return 'N'
    else:
        return 'Unknown'

def create_html_table(results, output_file):
    """
    Create HTML formatted frequency table
    """
    total_sequences = results['total_sequences']
    major_hg_counts = results['major_hg_counts']
    l_total = results['l_total']
    m_total = results['m_total'] 
    n_total = results['n_total']
    lineage_counts = results['lineage_counts']
    
    # Define haplogroup orders
    l_order = ['L3', 'L0', 'L2', 'L1', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9']
    m_order = ['M', 'D', 'C', 'G', 'E', 'Q', 'Z']
    n_order = ['H', 'U', 'B', 'J', 'T', 'K', 'F', 'A', 'R', 'HV', 'N', 'I', 'V', 'X', 'W', 'P', 'Y', 'S', 'O']
    
    def get_lineage_data(lineage_letter, haplogroup_order):
        """Get haplogroup data for a specific lineage"""
        lineage_data = []
        lineage_total = 0
        
        for hg in haplogroup_order:
            if hg in major_hg_counts.index:
                count = major_hg_counts[hg]
                lineage_total += count
                percentage = (count / lineage_total) * 100 if lineage_total > 0 else 0
                lineage_data.append((hg, count, percentage))
        
        # Add any haplogroups not in predefined order
        for hg in major_hg_counts.index:
            if categorize_haplogroup(hg) == lineage_letter and hg not in haplogroup_order:
                count = major_hg_counts[hg]
                lineage_total += count
                percentage = (count / lineage_total) * 100 if lineage_total > 0 else 0
                lineage_data.append((hg, count, percentage))
        
        # Recalculate percentages within lineage
        for i, (hg, count, _) in enumerate(lineage_data):
            percentage = (count / lineage_total) * 100 if lineage_total > 0 else 0
            lineage_data[i] = (hg, count, percentage)
            
        return lineage_data, lineage_total
    
    l_data, l_total = get_lineage_data('L', l_order)
    m_data, m_total = get_lineage_data('M', m_order)
    n_data, n_total = get_lineage_data('N', n_order)
    
    # Generate HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Haplogroup Frequency Analysis</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            text-align: center;
            color: #2c3e50;
            margin-bottom: 10px;
            font-size: 2.2em;
        }}
        .subtitle {{
            text-align: center;
            color: #7f8c8d;
            margin-bottom: 30px;
            font-size: 1.1em;
        }}
        .main-table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
        }}
        .lineage-header {{
            background: linear-gradient(135deg, #3498db, #2980b9);
            color: white;
            text-align: center;
            padding: 15px;
            font-size: 1.3em;
            font-weight: bold;
        }}
        .lineage-subheader {{
            background-color: #ecf0f1;
            text-align: center;
            padding: 10px;
            font-style: italic;
            color: #2c3e50;
        }}
        .column-header {{
            background-color: #34495e;
            color: white;
            padding: 10px;
            text-align: center;
            font-weight: bold;
            border-bottom: 2px solid #2c3e50;
        }}
        .data-cell {{
            padding: 8px 12px;
            text-align: center;
            border-bottom: 1px solid #bdc3c7;
        }}
        .data-cell:nth-child(1) {{ text-align: left; font-weight: bold; }}
        .data-cell:nth-child(2) {{ text-align: right; }}
        .data-cell:nth-child(3) {{ text-align: right; }}
        .total-row {{
            background-color: #f8f9fa;
            font-weight: bold;
            border-top: 2px solid #34495e;
        }}
        .overall-row {{
            background-color: #e8f4f8;
            font-weight: bold;
            color: #2c3e50;
        }}
        .lineage-column {{
            width: 33.33%;
            vertical-align: top;
            border-right: 2px solid #bdc3c7;
        }}
        .lineage-column:last-child {{
            border-right: none;
        }}
        .summary-section {{
            background-color: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            margin-top: 20px;
        }}
        .summary-title {{
            font-size: 1.4em;
            font-weight: bold;
            color: #2c3e50;
            margin-bottom: 15px;
        }}
        .summary-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }}
        .stat-box {{
            background-color: white;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid #3498db;
        }}
        .stat-value {{
            font-size: 1.8em;
            font-weight: bold;
            color: #2c3e50;
        }}
        .stat-label {{
            color: #7f8c8d;
            font-size: 0.9em;
            margin-top: 5px;
        }}
        .african {{ border-left-color: #e74c3c; }}
        .asian {{ border-left-color: #f39c12; }}
        .eurasian {{ border-left-color: #27ae60; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Haplogroup Frequency Analysis</h1>
        <p class="subtitle">Lineage distribution of {total_sequences:,} sequences from haplogrep3 analysis</p>
        
        <table class="main-table">
            <tr>
                <td class="lineage-column">
                    <table style="width: 100%;">
                        <tr><td class="lineage-header">L</td></tr>
                        <tr><td class="lineage-subheader">Lineages<br>"African"</td></tr>
                        <tr>
                            <td class="column-header" style="width: 25%;">hg</td>
                            <td class="column-header" style="width: 40%;">#</td>
                            <td class="column-header" style="width: 35%;">%</td>
                        </tr>"""
    
    # Add L lineage data
    for hg, count, percentage in l_data:
        html_content += f"""
                        <tr>
                            <td class="data-cell">{hg}</td>
                            <td class="data-cell">{count:,}</td>
                            <td class="data-cell">{percentage:.0f}%</td>
                        </tr>"""
    
    l_overall_pct = (l_total / total_sequences) * 100
    html_content += f"""
                        <tr class="total-row">
                            <td class="data-cell">Total</td>
                            <td class="data-cell">{l_total:,}</td>
                            <td class="data-cell">100%</td>
                        </tr>
                        <tr class="overall-row">
                            <td colspan="3" class="data-cell">Overall {l_overall_pct:.0f}% ({l_total:,} / {total_sequences:,})</td>
                        </tr>
                    </table>
                </td>
                
                <td class="lineage-column">
                    <table style="width: 100%;">
                        <tr><td class="lineage-header">M</td></tr>
                        <tr><td class="lineage-subheader">Lineages<br>"Asian"</td></tr>
                        <tr>
                            <td class="column-header" style="width: 25%;">hg</td>
                            <td class="column-header" style="width: 40%;">#</td>
                            <td class="column-header" style="width: 35%;">%</td>
                        </tr>"""
    
    # Add M lineage data
    for hg, count, percentage in m_data:
        html_content += f"""
                        <tr>
                            <td class="data-cell">{hg}</td>
                            <td class="data-cell">{count:,}</td>
                            <td class="data-cell">{percentage:.0f}%</td>
                        </tr>"""
    
    m_overall_pct = (m_total / total_sequences) * 100
    html_content += f"""
                        <tr class="total-row">
                            <td class="data-cell">Total</td>
                            <td class="data-cell">{m_total:,}</td>
                            <td class="data-cell">100%</td>
                        </tr>
                        <tr class="overall-row">
                            <td colspan="3" class="data-cell">Overall {m_overall_pct:.0f}% ({m_total:,} / {total_sequences:,})</td>
                        </tr>
                    </table>
                </td>
                
                <td class="lineage-column">
                    <table style="width: 100%;">
                        <tr><td class="lineage-header">N</td></tr>
                        <tr><td class="lineage-subheader">Lineages<br>"Eurasian"</td></tr>
                        <tr>
                            <td class="column-header" style="width: 25%;">hg</td>
                            <td class="column-header" style="width: 40%;">#</td>
                            <td class="column-header" style="width: 35%;">%</td>
                        </tr>"""
    
    # Add N lineage data
    for hg, count, percentage in n_data:
        html_content += f"""
                        <tr>
                            <td class="data-cell">{hg}</td>
                            <td class="data-cell">{count:,}</td>
                            <td class="data-cell">{percentage:.0f}%</td>
                        </tr>"""
    
    n_overall_pct = (n_total / total_sequences) * 100
    html_content += f"""
                        <tr class="total-row">
                            <td class="data-cell">Total</td>
                            <td class="data-cell">{n_total:,}</td>
                            <td class="data-cell">100%</td>
                        </tr>
                        <tr class="overall-row">
                            <td colspan="3" class="data-cell">Overall {n_overall_pct:.0f}% ({n_total:,} / {total_sequences:,})</td>
                        </tr>
                    </table>
                </td>
            </tr>
        </table>
        
        <div class="summary-section">
            <div class="summary-title">Summary Statistics</div>
            <div class="summary-stats">
                <div class="stat-box">
                    <div class="stat-value">{total_sequences:,}</div>
                    <div class="stat-label">Total Sequences Analyzed</div>
                </div>
                <div class="stat-box african">
                    <div class="stat-value">{l_total:,}</div>
                    <div class="stat-label">L lineages (African) - {l_overall_pct:.1f}%</div>
                </div>
                <div class="stat-box asian">
                    <div class="stat-value">{m_total:,}</div>
                    <div class="stat-label">M lineages (Asian) - {m_overall_pct:.1f}%</div>
                </div>
                <div class="stat-box eurasian">
                    <div class="stat-value">{n_total:,}</div>
                    <div class="stat-label">N lineages (Eurasian) - {n_overall_pct:.1f}%</div>
                </div>"""
    
    # Add unknown count if present
    unknown_count = lineage_counts.get('Unknown', 0)
    if unknown_count > 0:
        unknown_pct = (unknown_count / total_sequences) * 100
        html_content += f"""
                <div class="stat-box">
                    <div class="stat-value">{unknown_count:,}</div>
                    <div class="stat-label">Unknown/Unclassified - {unknown_pct:.1f}%</div>
                </div>"""
    
    html_content += f"""
            </div>
        </div>
        
        <div style="margin-top: 20px; text-align: center; color: #7f8c8d; font-size: 0.9em;">
            Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')} | 
            Unique major haplogroups: {len(major_hg_counts)}
        </div>
    </div>
</body>
</html>"""
    
    # Save HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"HTML frequency table saved to: {output_file}")

def create_frequency_table(csv_file_path, output_file=None, output_format='txt'):
    """
    Create frequency table from haplogrep3 CSV output
    
    Args:
        csv_file_path: Path to the CSV file
        output_file: Output file name (extension will be added based on format)
        output_format: 'txt', 'html', or 'both'
    """
    # Read the CSV file efficiently for large datasets
    try:
        # Read in chunks to handle large files efficiently
        chunk_size = 10000
        chunks = []
        
        print(f"Reading CSV file: {csv_file_path}")
        for chunk in pd.read_csv(csv_file_path, chunksize=chunk_size):
            chunks.append(chunk)
            
        df = pd.concat(chunks, ignore_index=True)
        print(f"Loaded {len(df)} records from {csv_file_path}")
        
        # Clean up memory
        del chunks
        
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    # Check if Haplogroup column exists
    if 'Haplogroup' not in df.columns:
        print("Error: 'Haplogroup' column not found in CSV file")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Extract major haplogroups
    df['MajorHaplogroup'] = df['Haplogroup'].apply(extract_major_haplogroup)
    df['Lineage'] = df['MajorHaplogroup'].apply(categorize_haplogroup)
    
    # Calculate frequencies
    total_sequences = len(df)
    lineage_counts = df['Lineage'].value_counts()
    major_hg_counts = df['MajorHaplogroup'].value_counts()
    
    # Prepare output
    output_lines = []
    output_lines.append(f"Lineage distribution of {total_sequences:,} sequences from haplogrep3 analysis:")
    output_lines.append("")
    
    # Define the order for each lineage
    l_order = ['L3', 'L0', 'L2', 'L1', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9']
    m_order = ['M', 'D', 'C', 'G', 'E', 'Q', 'Z']
    n_order = ['H', 'U', 'B', 'J', 'T', 'K', 'F', 'A', 'R', 'HV', 'N', 'I', 'V', 'X', 'W', 'P', 'Y', 'S', 'O']
    
    # Create three-column layout
    def create_lineage_section(lineage_name, lineage_desc, haplogroup_order, total_count):
        lines = []
        lines.append(f"        {lineage_name}")
        lines.append(f"     Lineages")
        lines.append(f'    "{lineage_desc}"')
        lines.append("")
        lines.append("  hg      #        %")
        
        lineage_total = 0
        for hg in haplogroup_order:
            if hg in major_hg_counts.index:
                count = major_hg_counts[hg]
                lineage_total += count
                percentage = (count / total_count) * 100
                lines.append(f"  {hg:<3}  {count:>6,}  {percentage:>4.0f}%")
        
        # Add any haplogroups not in the predefined order
        lineage_letter = lineage_name
        for hg in major_hg_counts.index:
            if categorize_haplogroup(hg) == lineage_letter and hg not in haplogroup_order:
                count = major_hg_counts[hg]
                lineage_total += count
                percentage = (count / total_count) * 100
                lines.append(f"  {hg:<3}  {count:>6,}  {percentage:>4.0f}%")
        
        lines.append("")
        lines.append(f" Total   {lineage_total:>6,}  100%")
        
        overall_percentage = (lineage_total / total_count) * 100
        lines.append(f"Overall {overall_percentage:.0f}% ({lineage_total:,} / {total_count:,})")
        
        return lines, lineage_total
    
    # Generate sections for each lineage
    l_lines, l_total = create_lineage_section('L', 'African', l_order, total_sequences)
    m_lines, m_total = create_lineage_section('M', 'Asian', m_order, total_sequences)
    n_lines, n_total = create_lineage_section('N', 'Eurasian', n_order, total_sequences)
    
    # Combine in three-column format (simplified for readability)
    max_lines = max(len(l_lines), len(m_lines), len(n_lines))
    
    # Pad shorter sections
    while len(l_lines) < max_lines:
        l_lines.append("")
    while len(m_lines) < max_lines:
        m_lines.append("")
    while len(n_lines) < max_lines:
        n_lines.append("")
    
    # Print side by side
    for i in range(max_lines):
        line = f"{l_lines[i]:<35} {m_lines[i]:<35} {n_lines[i]}"
        output_lines.append(line)
    
    # Add summary statistics
    output_lines.append("")
    output_lines.append("="*100)
    output_lines.append("SUMMARY STATISTICS:")
    output_lines.append(f"Total sequences analyzed: {total_sequences:,}")
    output_lines.append(f"L lineages (African): {l_total:,} ({(l_total/total_sequences)*100:.1f}%)")
    output_lines.append(f"M lineages (Asian): {m_total:,} ({(m_total/total_sequences)*100:.1f}%)")
    output_lines.append(f"N lineages (Eurasian): {n_total:,} ({(n_total/total_sequences)*100:.1f}%)")
    
    # Handle unknown/unclassified
    unknown_count = lineage_counts.get('Unknown', 0)
    if unknown_count > 0:
        output_lines.append(f"Unknown/Unclassified: {unknown_count:,} ({(unknown_count/total_sequences)*100:.1f}%)")
    
    # Display output
    for line in output_lines:
        print(line)
    
    # Prepare results dictionary
    results = {
        'total_sequences': total_sequences,
        'lineage_counts': lineage_counts,
        'major_hg_counts': major_hg_counts,
        'l_total': l_total,
        'm_total': m_total,
        'n_total': n_total,
        'dataframe': df
    }
    
    # Save to file(s) based on format
    if output_file:
        base_name = output_file.rsplit('.', 1)[0] if '.' in output_file else output_file
        
        if output_format in ['txt', 'both']:
            txt_file = f"{base_name}.txt"
            with open(txt_file, 'w') as f:
                f.write('\n'.join(output_lines))
            print(f"\nText frequency table saved to: {txt_file}")
        
        if output_format in ['html', 'both']:
            html_file = f"{base_name}.html"
            create_html_table(results, html_file)
    
    # Return detailed breakdown for further analysis
    return results

def analyze_haplogroup_diversity(results):
    """
    Additional analysis of haplogroup diversity
    """
    if not results:
        return
    
    df = results['dataframe']
    
    print("\n" + "="*50)
    print("DETAILED HAPLOGROUP ANALYSIS:")
    print("="*50)
    
    # Top 20 most common haplogroups
    print("\nTop 20 most common major haplogroups:")
    top_20 = results['major_hg_counts'].head(20)
    for i, (hg, count) in enumerate(top_20.items(), 1):
        percentage = (count / results['total_sequences']) * 100
        print(f"{i:2d}. {hg:<4} {count:>6,} ({percentage:>5.1f}%)")
    
    # Diversity metrics
    print(f"\nHaplogroup diversity metrics:")
    print(f"Total unique major haplogroups: {len(results['major_hg_counts'])}")
    
    # Calculate Shannon diversity index
    import numpy as np
    frequencies = results['major_hg_counts'].values / results['total_sequences']
    shannon_diversity = -np.sum(frequencies * np.log(frequencies))
    print(f"Shannon diversity index: {shannon_diversity:.3f}")
    
    # Simpson's diversity index
    simpson_diversity = 1 - np.sum(frequencies**2)
    print(f"Simpson's diversity index: {simpson_diversity:.3f}")

if __name__ == "__main__":
    import sys
    import argparse
    
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Generate haplogroup frequency table from haplogrep3 output')
    parser.add_argument('input_file', nargs='?', default='all_haplogroups.csv',
                       help='Input CSV file from haplogrep3 (default: all_haplogroups.csv)')
    parser.add_argument('-o', '--output', default='haplogroup_frequency_table',
                       help='Output file base name (default: haplogroup_frequency_table)')
    parser.add_argument('-f', '--format', choices=['txt', 'html', 'both'], default='txt',
                       help='Output format: txt, html, or both (default: txt)')
    
    # Parse arguments (fall back to simple sys.argv if argparse fails)
    try:
        args = parser.parse_args()
        csv_file_path = args.input_file
        output_file = args.output
        output_format = args.format
    except:
        # Fallback for simple usage
        csv_file_path = sys.argv[1] if len(sys.argv) > 1 else "all_haplogroups.csv"
        output_file = "haplogroup_frequency_table"
        output_format = "txt"
    
    # Create frequency table
    print("Processing haplogroup data...")
    print(f"Input file: {csv_file_path}")
    print(f"Output format: {output_format}")
    
    results = create_frequency_table(csv_file_path, output_file=output_file, output_format=output_format)
    
    if results:
        # Additional analysis
        analyze_haplogroup_diversity(results)
        
        print(f"\nProcessing complete!")
        print(f"Input file: {csv_file_path}")
        if output_format == 'both':
            print(f"Output files: {output_file}.txt, {output_file}.html")
        elif output_format == 'html':
            print(f"Output file: {output_file}.html")
        else:
            print(f"Output file: {output_file}.txt")
        
        if output_format in ['html', 'both']:
            print(f"\nOpen the HTML file in your web browser to view the formatted table!")
    else:
        print("Error: Failed to process the data.")
