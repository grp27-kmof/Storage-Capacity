import os
import csv
import re
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from pymatgen.analysis.bond_valence import BVAnalyzer


INPUT_DIR = "cif"
OUTPUT_CIF_DIR = "output_cif"
OUTPUT_CSV = "results.csv"
os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_CIF_DIR, exist_ok=True)

# Fundamental Constants
FARADAY = 96485.33  # C/mol (Charge per mole of electrons)
MIN_K_HOST_DIST = 2.00  # Minimum K-Host atom distance (Angstroms)
MIN_K_K_DIST = 3.50  # Updated to 3.50 Å to account for K+-K+ coulombic repulsion
GRAPHITE_K = 279.0  # Industry benchmark for KC8 (Graphite Anode)
K_VOL_EST = 11.0  # Estimated volume footprint of K+ in a lattice (A^3)

# Stoichiometric Mapping for Alloys and Conversion Materials
K_MECHANISM_MAP = {
    'P': 1.0, 'Sn': 1.0, 'Bi': 3.0, 'Sb': 3.0, 'Ge': 1.0, 
    'Si': 1.0, 'Pb': 1.0, 'S': 2.0, 'Se': 2.0, 'Te': 2.0,
    'Al': 1.0, 'Mg': 1.0, 'Ca': 1.0, 'Zn': 1.0, 'Sr': 1.0, 'Ba': 1.0
}

def load_structure_robust(cif_path):
    """
    Loads a CIF, tolerating '?' placeholders in numeric fields
    (unrefined occupancy / U_iso / coordinate values), which crash
    the default pymatgen float parser.
    """
    try:
        return Structure.from_file(cif_path)
    except ValueError as e:
        if "'?'" not in str(e) and "convert string to float: '?'" not in str(e):
            raise
        print(f"  [PARSE WARN] Malformed numeric field ('?') detected -- sanitizing and retrying.")
        with open(cif_path, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        # Replace standalone '?' tokens (whitespace-delimited) with 0.0
        content_clean = re.sub(r'(?<=\s)\?(?=\s|\n|$)', '0.0', content)
        parser = CifParser.from_str(content_clean)
        structures = parser.parse_structures(primitive=False)
        return structures[0]

def find_k_insertion_sites(calc, host_struct, n_final, energy_threshold=0.2):
    """
    Greedily samples stable BVSE voxels across open paths. This prevents
    the center-of-mass bug by allowing multiple distinct ions to be placed
    sequentially along a single continuous low-energy channel.
    """
    grid = calc.grid_data
    shape = np.array(grid.shape)
    lattice = host_struct.lattice
    host_fracs = np.array([s.frac_coords for s in host_struct])

    # Dynamic threshold relaxation loop
    threshold = energy_threshold
    for _ in range(5):
        stable_indices = np.argwhere(grid < threshold)
        if len(stable_indices) > 0:
            break
        threshold *= 2.0
        print(f"  [WARN] No sites below threshold; relaxing to {threshold:.3f}")

    if len(stable_indices) == 0:
        return []

    # Sort candidates by energy level (lowest energy valleys filled first)
    energies = grid[tuple(stable_indices.T)]
    stable_indices = stable_indices[np.argsort(energies)]

    accepted = []
    for idx in stable_indices:
        frac = idx / shape
        
        # Periodic Boundary Condition (PBC) distance check against Host framework
        diff_host = frac - host_fracs
        diff_host -= np.round(diff_host)
        cart_host = diff_host @ lattice.matrix
        if np.any(np.linalg.norm(cart_host, axis=1) < MIN_K_HOST_DIST):
            continue

        # PBC distance check against previously accepted Potassium ions
        if accepted:
            diff_k = frac - np.array(accepted)
            diff_k -= np.round(diff_k)
            cart_k = diff_k @ lattice.matrix
            if np.any(np.linalg.norm(cart_k, axis=1) < MIN_K_K_DIST):
                continue

        accepted.append(frac)
        if len(accepted) == n_final:
            break

    return accepted

def find_k_sites_geometric(host_struct, n_final, resolution=0.5):
    lattice = host_struct.lattice
    host_fracs = np.array([s.frac_coords for s in host_struct])

    min_host_dist = MIN_K_HOST_DIST
    res = resolution
    accepted = []

    for attempt in range(6):
        na = max(3, int(lattice.a / res))
        nb = max(3, int(lattice.b / res))
        nc = max(3, int(lattice.c / res))

        va = np.linspace(0, 1, na, endpoint=False)
        vb = np.linspace(0, 1, nb, endpoint=False)
        vc = np.linspace(0, 1, nc, endpoint=False)

        grid_a, grid_b, grid_c = np.meshgrid(va, vb, vc, indexing='ij')
        grid_fracs = np.stack([grid_a.flatten(), grid_b.flatten(), grid_c.flatten()], axis=-1)

        diff = grid_fracs[:, None, :] - host_fracs[None, :, :]
        diff -= np.round(diff)  # Apply Periodic Boundary Conditions
        cart_diff = diff @ lattice.matrix
        dists = np.linalg.norm(cart_diff, axis=2)
        min_dists_to_host = np.min(dists, axis=1)

        valid_indices = np.where(min_dists_to_host >= min_host_dist)[0]

        if len(valid_indices) > 0:
            void_candidates = list(zip(min_dists_to_host[valid_indices], grid_fracs[valid_indices]))
            void_candidates.sort(key=lambda x: -x[0])  # Prioritize widest open voids

            accepted = []
            for _, frac in void_candidates:
                if accepted:
                    diff_k = frac - np.array(accepted)
                    diff_k -= np.round(diff_k)
                    if np.any(np.linalg.norm(diff_k @ lattice.matrix, axis=1) < MIN_K_K_DIST):
                        continue
                accepted.append(frac)
                if len(accepted) == n_final:
                    break

            if accepted:
                if attempt > 0:
                    print(f"  [SITE INFO] Dense host: relaxed min_host_dist={min_host_dist:.2f} A, "
                          f"resolution={res:.2f} A to find {len(accepted)} site(s).")
                return accepted

        # Nothing valid at this tolerance -- tighten the grid and relax the
        # host-clash distance floor (never below 1.0 A, which is already
        # unphysically close and signals a genuinely non-porous host).
        res = max(0.1, res * 0.5)
        min_host_dist = max(1.0, min_host_dist * 0.85)

    return accepted

def calculate_universal_k_capacity(cif_path):
    """
    Core Evaluation Engine. Resolves structural limits and scales 
    calculations directly against full unit cell mass boundaries.
    """
    print(f"--- Comprehensive Potassium Analysis: {os.path.basename(cif_path)} ---")
    
    struct = load_structure_robust(cif_path)
    host_struct = struct.copy()
    
    # Isolate empty framework matrix
    if "K" in [s.specie.symbol for s in struct if s.is_ordered]:
        host_struct.remove_species(["K"])
    
    mw_host = host_struct.composition.weight
    formula = host_struct.composition.reduced_formula
    vol_cell = host_struct.volume
    
    # Check framework classification
    has_carbon = any(site.is_ordered and site.specie.symbol == "C" for site in host_struct)

    # -------------------------------------------------------------------------
    # STAGE 2: REDOX RESOLUTION (n_redox)
    # -------------------------------------------------------------------------
    n_redox = 0
    n_zr = sum(1 for s in host_struct if s.is_ordered and s.specie.symbol == "Zr")

    if n_zr >= 6 and has_carbon:
        # Benchmark for UiO-66 framework topologies
        n_redox = (n_zr / 6.0) * 40.0
    else:
        composition_dict = host_struct.composition.as_dict()
        for site in host_struct:
            if not site.is_ordered:
                continue

            sym = site.specie.symbol
            if sym in K_MECHANISM_MAP:
                n_redox += K_MECHANISM_MAP[sym]
            elif site.specie.is_transition_metal:
                n_redox += 1.0  # Normalized to 1e- standard framework transition reduction
        
        # Corrected Organic Linker Redox allocation
        if has_carbon:
            num_oxygen = composition_dict.get("O", 0)
            n_redox += (num_oxygen / 2.0)  # Count redox pairs (carboxyls/carbonyls)

    # -------------------------------------------------------------------------
    # STAGE 3: GEOMETRIC SPACE SCREENING (n_geo)
    # -------------------------------------------------------------------------
    bvse_success = False
    try:
        from bvlain import Lain
        calc = Lain(verbose=False)
        calc.read_structure(host_struct)
        params = {'mobile_ion': 'K1+', 'r_cut': 10.0, 'resolution': 0.2, 'k': 100}
        calc.bvse_distribution(**params)
        
        stable_voxels = np.sum(calc.grid_data < 0.2)
        n_geo = (stable_voxels * (0.2**3)) / K_VOL_EST
        bvse_success = True
    except Exception as e:
        print(f"  [BVSE WARN] {e} -- routing to geometric mesh grid finder.")
        void_fraction = 0.5 if has_carbon else 0.1
        n_geo = (vol_cell * void_fraction) / K_VOL_EST

    # -------------------------------------------------------------------------
    # STAGE 4: INTERSECTION ANALYSIS & MAPPING
    # -------------------------------------------------------------------------
    n_final = min(n_redox, n_geo) if has_carbon else n_redox
    n_final = max(1, int(round(n_final)))

    # NOTE: Site placement is now attempted for EVERY host (carbon or not).
    # Previously this was gated behind `if has_carbon:`, which silently
    # skipped inorganic/alloy frameworks (e.g. CoP) and any framework
    # where pymatgen didn't resolve explicit carbon sites (e.g. simplified
    # Prussian Blue CIFs).
    if bvse_success:
        insertion_sites = find_k_insertion_sites(calc, host_struct, n_final)
    else:
        insertion_sites = find_k_sites_geometric(host_struct, n_final)

    if not insertion_sites:
        print(f"  [SITE WARN] 0 sites found for n_final={n_final} "
              f"(host elements: {sorted(set(s.specie.symbol for s in host_struct if s.is_ordered))})")

    if insertion_sites:
        base_name = os.path.splitext(os.path.basename(cif_path))[0]
        inserted_struct = host_struct.copy()
        for site in insertion_sites:
            inserted_struct.append("K", site, validate_proximity=False)
        output_cif = os.path.join(OUTPUT_CIF_DIR, f"{base_name}.cif")
        inserted_struct.to(filename=output_cif)
        print(f"  K-inserted CIF successfully generated: {output_cif}")

    # Faraday Capacity Conversion using absolute Unit Cell Mass
    capacity = (n_final * FARADAY) / (3.6 * mw_host)
    improvement = ((capacity / GRAPHITE_K) - 1) * 100

    return {
        "cif_file": os.path.basename(cif_path),
        "Formula": formula,
        "Molar_Mass_Host": round(mw_host, 2),
        "Max_K_Ions_n": int(n_final),
        "Theoretical_Capacity_mAh_g": round(capacity, 2),
        "Improvement_vs_KGr_perc": round(improvement, 2)
    }

def save_results(data, filename):
    file_exists = os.path.isfile(filename)
    with open(filename, 'a', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=data.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(data)

if __name__ == "__main__":
    cif_files = sorted(
        os.path.join(INPUT_DIR, f) for f in os.listdir(INPUT_DIR) if f.lower().endswith(".cif")
    )

    if not cif_files:
        print(f"No target files isolated inside '{INPUT_DIR}' directory.")
    else:
        if os.path.isfile(OUTPUT_CSV):
            os.remove(OUTPUT_CSV)
            
        for cif_path in cif_files:
            try:
                res = calculate_universal_k_capacity(cif_path)
                save_results(res, OUTPUT_CSV)
            except Exception as e:
                print(f"[ERROR] Failed processing targets on node: {cif_path}\n{e}")