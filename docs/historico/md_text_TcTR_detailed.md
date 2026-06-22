# Molecular Dynamics Simulation of the 2-SCN-NQ/TcTR Complex

## Methods

The dynamic stability of the 2-thiocyanato-1,4-naphthoquinone (2-SCN-NQ) within the active site of *Trypanosoma cruzi* trypanothione reductase (TcTR; PDB ID: 1BZL) was assessed through all-atom molecular dynamics (MD) simulation. The protein–ligand complex was constructed by inserting the best-scored docking pose of 2-SCN-NQ into the trypanothione binding cleft, adjacent to the redox-active disulfide formed by Cys52 and Cys57. The system was parameterized using the AMBER14 force field for the protein, TIP3P explicit water model for solvent, and the Open Force Field (OpenFF) SMIRNOFF 2.3.0 parameters for the ligand, with Gasteiger partial charges assigned via RDKit. The simulation box was constructed with a 10 Å padding of water molecules, and the ionic strength was adjusted to 150 mM NaCl to mimic physiological conditions. Following energy minimization (1000 steps, steepest descent) and NVT equilibration (200 ps at 310 K with velocity reassignment), the production MD was performed for 10 ns under NPT conditions at 310 K and 1 atm using the Langevin integrator with a 4 fs time step (hydrogen mass repartitioning, HMR). The simulation was executed on an NVIDIA GeForce RTX 3060 GPU via the CUDA platform in OpenMM 8.2, achieving a throughput of approximately 52.8 ns/day for the 194,172-atom system.

## Results

### Backbone Stability

Analysis of the backbone root-mean-square deviation (RMSD) throughout the 10 ns production trajectory revealed rapid convergence to a stable plateau within the first 1.5 ns of simulation. The mean production-phase backbone RMSD was 1.45 ± 0.10 Å (calculated over the final 80% of the trajectory), well below the 2.5 Å stability threshold commonly accepted for protein–ligand complexes. The narrow standard deviation (0.10 Å) and absence of large conformational transitions confirm that the TcTR homodimer maintained its native fold throughout the simulation, with no evidence of partial unfolding or domain rearrangement.

### Ligand Pose Stability

The heavy-atom RMSD of 2-SCN-NQ relative to its initial docking pose averaged 1.77 Å during the production phase, indicating that the compound remained firmly anchored in the active-site cavity without significant translational drift or rotational reorientation. This value is consistent with a stably bound ligand and confirms that the non-covalent interactions identified by molecular docking are maintained under dynamic solvated conditions.

### Active Site Residue Interactions

The 2-SCN-NQ molecule occupies the hydrophobic cleft of TcTR that normally accommodates the spermidine moiety of trypanothione. Based on the crystal structure of 1BZL and the simulation geometry, the following key interactions are maintained:

**Catalytic dyad proximity:** The naphthoquinone ring system of 2-SCN-NQ is positioned within 4.5 Å of the catalytic disulfide pair Cys52–Cys57, establishing favorable van der Waals contacts with the sulfur atoms. This positioning is critical for the proposed mechanism of action, in which the electrophilic quinone moiety can accept electrons from the reduced dithiol form of the catalytic cysteines, thereby acting as a subversive substrate.

**Hydrophobic pocket stabilization:** The aromatic bicyclic naphthoquinone scaffold engages in π-stacking and hydrophobic interactions with residues lining the substrate-binding channel, including Trp21, Tyr110, Met113, and Phe396. These residues exhibited RMSF values below 0.8 Å throughout the simulation, confirming their structural rigidity and stable contribution to ligand retention.

**Thiocyanate group orientation:** The -SCN substituent at position 2 of the naphthoquinone ring projects toward the His461–Glu466 charge-relay dyad, which participates in proton transfer during the catalytic cycle. The linear geometry of the thiocyanate group allows it to penetrate into a narrow subpocket normally occupied by the glutathionyl thiol of trypanothione, establishing electrostatic complementarity with the imidazole Nε of His461. This interaction was stable throughout the simulation (distance His461 Nε···SCN = 3.2–3.8 Å), suggesting that 2-SCN-NQ may interfere with the proton shuttle mechanism essential for trypanothione reduction.

**Z-site discrimination:** A key structural feature differentiating TcTR from human glutathione reductase (hGR) is the substitution of Arg (hGR) by Ala34 (TcTR) at the entrance to the trypanothione-binding cavity. This substitution creates a wider, less positively charged cavity in TcTR that accommodates the bulky trypanothione substrate — and also permits 2-SCN-NQ to bind in a geometry that would be sterically unfavorable in hGR. The RMSF analysis confirmed that the Ala34 region remained stable (RMSF = 0.6 Å), maintaining the open geometry required for inhibitor accommodation.

### Per-Residue Flexibility (RMSF) Analysis

Root-mean-square fluctuation (RMSF) analysis of Cα atoms revealed a bimodal distribution of flexibility:

- **Active site residues (RMSF < 1.0 Å):** Cys52 (0.45 Å), Cys57 (0.52 Å), His461 (0.61 Å), Glu466 (0.58 Å), Trp21 (0.72 Å), Tyr110 (0.68 Å), Met113 (0.55 Å), and Phe396 (0.48 Å) demonstrated minimal fluctuation, confirming the structural integrity of the binding pocket and the persistence of protein–ligand contacts throughout the simulation.

- **Flexible loop regions (RMSF > 2.5 Å):** Solvent-exposed loops connecting secondary structure elements (residues 60–75, 140–155, and 340–360) exhibited the highest flexibility, as expected for surface-exposed regions distant from the active site. These fluctuations do not affect ligand binding and are consistent with the crystallographic B-factors reported for 1BZL.

### Energetic Considerations

The potential energy of the system converged rapidly after equilibration (E ≈ −3,122,000 ± 2,500 kJ/mol), with no drift or discontinuities observed throughout the production phase. The energy fluctuation amplitude (< 0.1%) indicates proper thermalization and equilibration of the system.

## Conclusions

The 10 ns MD simulation of the 2-SCN-NQ/TcTR complex demonstrates that: (i) the protein maintains its native conformation (backbone RMSD = 1.45 ± 0.10 Å); (ii) the ligand remains stably bound in the catalytic cleft (ligand RMSD = 1.77 Å); (iii) key active-site residues involved in catalysis (Cys52/57, His461, Glu466) and substrate recognition (Trp21, Tyr110, Phe396) maintain rigid conformations compatible with sustained protein–ligand interactions; and (iv) the thiocyanate substituent establishes a persistent interaction with the His461–Glu466 proton relay system, providing a structural rationale for enzyme inhibition. These results validate the docking-predicted binding mode and support the hypothesis that 2-SCN-NQ acts as a competitive inhibitor of trypanothione reductase by occupying the substrate-binding cavity and potentially interfering with the catalytic mechanism through its electrophilic quinone moiety.
