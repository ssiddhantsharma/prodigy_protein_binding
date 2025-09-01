#!/usr/bin/env python3

import os
import subprocess
import csv
import re
from multiprocessing import Pool
import argparse
from typing import Dict, List, Tuple, Optional

# Define binding strength thresholds
STRONG_BINDING = -9.0    # kcal/mol - Strong/good binding
MODERATE_BINDING = -7.0  # kcal/mol - Fair binding
SYMMETRIC_THRESHOLD = 1.5 # Max difference between A-B and A-C binding energy

def run_prodigy(pdb_path: str, chain1: str, chain2: str, temp: float = 25.0) -> Tuple[Optional[float], str, Dict]:
    """Run PRODIGY and extract all metrics."""
    basename = os.path.basename(pdb_path)
    interface = f"{chain1}-{chain2}"
    print(f"Analyzing {basename} ({interface})...")
    
    cmd = [
        '/nfs/scistore20/praetgrp/ssharma/.local/bin/prodigy',
        pdb_path,
        '--selection', chain1, chain2,
        '--temperature', str(temp)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check for "No contacts" error
        if result.returncode != 0:
            if "No contacts found for selection" in result.stderr:
                print(f"  {basename} ({interface}): No contacts")
                return None, 'NO_CONTACTS', {}
            else:
                print(f"  Error: {result.stderr[:100]}...")
                return None, 'ERROR', {}
        
        # Parse all metrics from the PRODIGY output
        metrics = {}
        
        # Extract binding affinity (ΔG) in kcal/mol
        dg_match = re.search(r'Predicted binding affinity \(kcal\.mol-1\):\s+([-\d.]+)', result.stdout)
        if dg_match:
            binding_energy = float(dg_match.group(1))
            metrics['binding_energy'] = binding_energy
            print(f"  {basename} ({interface}): {binding_energy:.2f} kcal/mol")
        else:
            print(f"  Could not extract binding energy from output")
            return None, 'PARSE_ERROR', {}
            
        # Extract dissociation constant (Kd) in M
        kd_match = re.search(r'Predicted dissociation constant \(M\) at .+?˚C:\s+([\d.e-]+)', result.stdout)
        if kd_match:
            metrics['kd'] = kd_match.group(1)
            
        # Extract total intermolecular contacts
        contacts_match = re.search(r'No\. of intermolecular contacts: (\d+)', result.stdout)
        if contacts_match:
            metrics['total_contacts'] = int(contacts_match.group(1))
            
        # Extract specific contact types
        contact_patterns = {
            'charged_charged': r'No\. of charged-charged contacts: (\d+)',
            'charged_polar': r'No\. of charged-polar contacts: (\d+)',
            'charged_apolar': r'No\. of charged-apolar contacts: (\d+)',
            'polar_polar': r'No\. of polar-polar contacts: (\d+)',
            'polar_apolar': r'No\. of apolar-polar contacts: (\d+)',
            'apolar_apolar': r'No\. of apolar-apolar contacts: (\d+)'
        }
        
        for contact_type, pattern in contact_patterns.items():
            match = re.search(pattern, result.stdout)
            if match:
                metrics[contact_type] = int(match.group(1))
        
        # Extract NIS percentages
        apolar_nis_match = re.search(r'Percentage of apolar NIS residues: ([\d.]+)', result.stdout)
        if apolar_nis_match:
            metrics['perc_apolar_nis'] = float(apolar_nis_match.group(1))
            
        charged_nis_match = re.search(r'Percentage of charged NIS residues: ([\d.]+)', result.stdout)
        if charged_nis_match:
            metrics['perc_charged_nis'] = float(charged_nis_match.group(1))
        
        return binding_energy, 'SUCCESS', metrics
            
    except Exception as e:
        print(f"  Exception: {str(e)[:100]}...")
        return None, 'EXCEPTION', {}

def assess_structure_quality(ab_energy: float, ac_energy: float, energy_diff: float, bc_detected: bool) -> Tuple[str, str]:
    """Assess structure quality based on binding energies and symmetry."""
    quality_issues = []
    
    # Check binding strength for A-B interface
    if ab_energy > MODERATE_BINDING:
        quality_issues.append('weak_AB_binding')
    elif ab_energy > STRONG_BINDING:
        quality_issues.append('moderate_AB_binding')
        
    # Check binding strength for A-C interface
    if ac_energy > MODERATE_BINDING:
        quality_issues.append('weak_AC_binding')
    elif ac_energy > STRONG_BINDING:
        quality_issues.append('moderate_AC_binding')
        
    # Check symmetry
    if energy_diff > SYMMETRIC_THRESHOLD:
        quality_issues.append('asymmetric_binding')
        
    # Flag BC binding as a critical issue
    if bc_detected:
        quality_issues.append('BC_binding_detected')
        
    # Overall quality assessment
    if not quality_issues:
        quality_assessment = 'Good'
    elif 'BC_binding_detected' in quality_issues:
        quality_assessment = 'Poor (B-C binding)'
    elif 'weak_AB_binding' in quality_issues or 'weak_AC_binding' in quality_issues:
        quality_assessment = 'Poor (weak binding)'
    elif 'asymmetric_binding' in quality_issues:
        quality_assessment = 'Fair (asymmetric binding)'
    elif 'moderate_AB_binding' in quality_issues or 'moderate_AC_binding' in quality_issues:
        quality_assessment = 'Fair (moderate binding)'
    else:
        quality_assessment = 'Unknown'
        
    quality_issues_str = ','.join(quality_issues) if quality_issues else 'None'
    
    return quality_assessment, quality_issues_str

def process_structure(args: Tuple[str, float]) -> Dict:
    """Process a single PDB structure."""
    pdb_path, temp = args
    basename = os.path.basename(pdb_path)
    
    try:
        # Run PRODIGY for each interface
        ab_energy, ab_status, ab_metrics = run_prodigy(pdb_path, 'A', 'B', temp)
        ac_energy, ac_status, ac_metrics = run_prodigy(pdb_path, 'A', 'C', temp)
        bc_energy, bc_status, bc_metrics = run_prodigy(pdb_path, 'B', 'C', temp)
        
        # Create result dictionary
        result = {'pdb': basename}
        
        # Add interface energies
        result['AB_energy'] = ab_energy
        result['AC_energy'] = ac_energy
        
        # Add AB metrics
        for key, value in ab_metrics.items():
            result[f'AB_{key}'] = value
            
        # Add AC metrics
        for key, value in ac_metrics.items():
            result[f'AC_{key}'] = value
        
        # Handle BC binding detection
        bc_binding_detected = bc_status == 'SUCCESS' and bc_energy is not None
        result['BC_binding_detected'] = bc_binding_detected
        if bc_binding_detected:
            result['BC_energy'] = bc_energy
        
        # Calculate derived metrics for valid AB-AC pairs
        if ab_energy is not None and ac_energy is not None:
            energy_diff = abs(ab_energy - ac_energy)
            result['energy_diff'] = energy_diff
            result['avg_binding_energy'] = (ab_energy + ac_energy) / 2
            
            # Assess quality
            quality_assessment, quality_issues = assess_structure_quality(
                ab_energy, ac_energy, energy_diff, bc_binding_detected
            )
            result['quality_assessment'] = quality_assessment
            result['quality_issues'] = quality_issues
        
        return result
        
    except Exception as e:
        print(f"Error processing {basename}: {str(e)}")
        return {
            'pdb': basename,
            'error': str(e),
            'status': 'PROCESSING_ERROR'
        }

def categorize_and_sort_results(results: List[Dict]) -> Tuple[List[Dict], Dict[str, List[Dict]]]:
    """Categorize and sort results by quality and binding strength."""
    # Filter valid results
    valid_results = [r for r in results if 'avg_binding_energy' in r and r['avg_binding_energy'] is not None]
    invalid_results = [r for r in results if 'avg_binding_energy' not in r or r['avg_binding_energy'] is None]
    
    # Initialize categories
    categories = {
        'good': [],
        'fair_moderate': [],
        'fair_asymmetric': [],
        'poor': []
    }
    
    if not valid_results:
        return invalid_results, categories
    
    # Categorize by quality
    categories = {
        'good': [r for r in valid_results if r.get('quality_assessment') == 'Good'],
        'fair_moderate': [r for r in valid_results if r.get('quality_assessment') == 'Fair (moderate binding)'],
        'fair_asymmetric': [r for r in valid_results if r.get('quality_assessment') == 'Fair (asymmetric binding)'],
        'poor': [r for r in valid_results if r.get('quality_assessment', '').startswith('Poor')]
    }
    
    # Sort each category by binding strength (strongest first)
    for category in categories.values():
        category.sort(key=lambda x: x['avg_binding_energy'])
    
    # Combine in order of quality
    sorted_results = (categories['good'] + categories['fair_moderate'] + 
                     categories['fair_asymmetric'] + categories['poor'] + invalid_results)
    
    return sorted_results, categories

def write_results_csv(sorted_results: List[Dict], output_file: str) -> None:
    """Write results to CSV file with proper headers."""
    if not sorted_results:
        print("No results to write to file")
        return
        
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Define column order and headers
    field_headers = {
        'pdb': 'pdb',
        'quality_assessment': 'quality_assessment',
        'quality_issues': 'quality_issues',
        'avg_binding_energy': 'avg_binding_energy (kcal/mol)',
        'energy_diff': 'energy_diff (kcal/mol)',
        'AB_energy': 'AB_energy (kcal/mol)',
        'AB_kd': 'AB_kd (M)',
        'AB_total_contacts': 'AB_total_contacts',
        'AB_charged_charged': 'AB_charged_charged',
        'AB_charged_polar': 'AB_charged_polar',
        'AB_charged_apolar': 'AB_charged_apolar',
        'AB_polar_polar': 'AB_polar_polar',
        'AB_polar_apolar': 'AB_polar_apolar',
        'AB_apolar_apolar': 'AB_apolar_apolar',
        'AB_perc_apolar_nis': 'AB_perc_apolar_nis (%)',
        'AB_perc_charged_nis': 'AB_perc_charged_nis (%)',
        'AC_energy': 'AC_energy (kcal/mol)',
        'AC_kd': 'AC_kd (M)',
        'AC_total_contacts': 'AC_total_contacts',
        'AC_charged_charged': 'AC_charged_charged',
        'AC_charged_polar': 'AC_charged_polar',
        'AC_charged_apolar': 'AC_charged_apolar',
        'AC_polar_polar': 'AC_polar_polar',
        'AC_polar_apolar': 'AC_polar_apolar',
        'AC_apolar_apolar': 'AC_apolar_apolar',
        'AC_perc_apolar_nis': 'AC_perc_apolar_nis (%)',
        'AC_perc_charged_nis': 'AC_perc_charged_nis (%)',
        'BC_binding_detected': 'BC_binding_detected',
        'BC_energy': 'BC_energy (kcal/mol)'
    }
    
    ordered_fields = list(field_headers.keys())
    headers = [field_headers[field] for field in ordered_fields]
    
    # Write CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        
        for result in sorted_results:
            row = [result.get(field, '') for field in ordered_fields]
            writer.writerow(row)
    
    print(f"\nResults saved to: {output_file}")

def print_summary(results: List[Dict], categories: Dict[str, List[Dict]]) -> None:
    """Print summary statistics and top structures."""
    bc_contacts = sum(1 for r in results if r.get('BC_binding_detected', False))
    
    print(f"\nSummary:")
    print(f"  Total structures analyzed: {len(results)}")
    print(f"  Good quality structures: {len(categories['good'])}")
    print(f"  Fair quality structures: {len(categories['fair_moderate']) + len(categories['fair_asymmetric'])}")
    print(f"  Poor quality structures: {len(categories['poor'])}")
    
    if bc_contacts > 0:
        print(f"  WARNING: {bc_contacts} structures have B-C binding (flagged in CSV)")
    else:
        print(f"  No B-C binding detected in any structures (ideal)")
    
    # Print top structures
    if categories['good']:
        print("\nTop good structures (recommended for experimental testing):")
        for i, result in enumerate(categories['good'][:5]):
            print(f"{i+1}. {result['pdb']}")
            print(f"   A-B: {result['AB_energy']:.2f} kcal/mol (Kd: {result.get('AB_kd', 'N/A')})")
            print(f"   A-C: {result['AC_energy']:.2f} kcal/mol (Kd: {result.get('AC_kd', 'N/A')})")
            print(f"   Energy diff: {result['energy_diff']:.2f} kcal/mol")
            print(f"   Avg binding energy: {result['avg_binding_energy']:.2f} kcal/mol")
            print(f"   A-B contacts: {result.get('AB_total_contacts', 'N/A')}")
            print(f"   A-C contacts: {result.get('AC_total_contacts', 'N/A')}")
            print("")

def batch_process_structures(input_dir: str, output_file: str, temp: float = 25.0, num_processes: int = 4) -> None:
    """Process all structures in input directory."""
    # Find PDB files
    pdb_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) 
                 if f.endswith(('.pdb', '.cif'))]
    
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return
    
    print(f"Found {len(pdb_files)} structure files to analyze")
    print(f"Using binding thresholds:")
    print(f"  Good binding: below {STRONG_BINDING} kcal/mol")
    print(f"  Fair binding: between {STRONG_BINDING} and {MODERATE_BINDING} kcal/mol")
    print(f"  Symmetric binding threshold: {SYMMETRIC_THRESHOLD} kcal/mol maximum difference")
    
    # Process structures in parallel
    args_list = [(pdb, temp) for pdb in pdb_files]
    
    with Pool(processes=num_processes) as pool:
        results = pool.map(process_structure, args_list)
    
    # Filter out None results
    results = [r for r in results if r]
    
    # Sort and categorize results
    sorted_results, categories = categorize_and_sort_results(results)
    
    # Write CSV output
    write_results_csv(sorted_results, output_file)
    
    # Print summary
    print_summary(results, categories)

def main():
    """Main function to parse arguments and run analysis."""
    parser = argparse.ArgumentParser(description='Analyze binding affinities in protein complexes using PRODIGY.')
    parser.add_argument('--input', required=True, help='Directory containing PDB files')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--temp', type=float, default=25.0, help='Temperature in Celsius (default: 25.0)')
    parser.add_argument('--processes', type=int, default=4, help='Number of parallel processes (default: 4)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not os.path.isdir(args.input):
        print(f"Error: Input directory '{args.input}' does not exist")
        return
    
    if not args.output.endswith('.csv'):
        print("Warning: Output file should have .csv extension")
    
    # Run analysis
    batch_process_structures(args.input, args.output, args.temp, args.processes)

if __name__ == "__main__":
    main()
