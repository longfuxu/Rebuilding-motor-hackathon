# AF3 ligand-state control (apo / ADP / ATP, template OFF) — clean predictor-limit result

User ran native ring (5 gp16 + 5 Mg) with apo / +ADP / +ATP, **all with useStructureTemplate=false**. Ring geometry:

| state | R146 (per subunit) | radius | planarity_rms | max out-of-plane | ring shape |
|---|---|---|---|---|---|
| apo + Mg | ~7.2 | 27.2 Å | 0.08 Å | 0.1 Å | closed planar |
| ADP + Mg | ~8.6 | 28.0 Å | 0.07 Å | 0.1 Å | closed planar |
| **ATP + Mg (no template)** | ~7.6 | 27.9 Å | **0.01 Å** | **0.0 Å** | **closed planar** |

**All three are flat closed rings — even ATP with the template OFF does NOT open the ring.** So the template was
not the cause; **AF3 fundamentally does not model the ATP-driven allosteric opening**. This is now a rigorous
3-state control (apo/ADP/ATP): the predictor gives the same closed planar ring regardless of nucleotide. Clean,
honest, and it *justifies* the NMA / occupancy-ENM / Gō-MD line (the right tools for the conformational question).
7JQQ (experimental, 3-ATP) remains the ground truth for the open state. Structures: `structures/af3_sweep/ligand_states/`.
