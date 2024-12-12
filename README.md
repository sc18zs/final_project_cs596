# CSCI596 Final Project: Methods to Enhance the Parallel Molecular Dynamics (MD) Program

## Group Members
- Zihan Shen
- Yuzhu Tang
- Sichen Zhang
- Zhengtao Yao

## Problem Statement
Molecular Dynamics (MD) simulations are computationally intensive, requiring efficient strategies to address challenges related to memory access, parallel communication, and processing power. This project aims to identify and implement methods to enhance the performance of MD simulations.

Efficient cache utilization is paramount for achieving high computational speeds and energy efficiency in applications ranging from scientific simulations to machine learning. Poor cache performance can lead to significant bottlenecks, limiting the scalability and responsiveness of systems.

### Key Challenges
1. **Inefficient Memory Access Due to Triple Nested Loops**
   - The current implementation uses a triple loop to traverse each cell and accesses neighboring cells via nested loops to compute particle interactions.
   - When the processor processes cells, memory access is often non-contiguous, and the order of accessing different cells may span large memory regions, increasing cache miss rates.

3. **Communication Overhead in Parallel Computing**  
   - MD simulations divide the computational domain among multiple processors, requiring frequent information exchange between adjacent regions.
   - Non-optimized cell access order increases cross-node communication frequency, significantly impacting computation efficiency.

4. **Load Imbalance**  
   - Uneven particle distribution across processors leads to an imbalanced computational load.
   - Some processors handle significantly more calculations than others, reducing resource utilization and overall performance.

## Objectives
1. Use Hilbert and Morton curves to improve cache hit rates, reduce memory access overhead, and enhance computational efficiency. 
> + Improve memory access locality and minimize cache misses by reordering particles based on spatial locality.
> + Use performance analysis tools to monitor cache hit rates and overall execution time. Compare results with the baseline (non-optimized approach) to verify the impact of reordering.
> + Develop algorithms that map multi-dimensional data to one-dimensional memory addresses using Hilbert and Morton curves.

2. Using load balancing to ensure the computational work is evenly distributed across all processors during a molecular dynamics simulation. When some processors have more atoms to compute than others, it creates inefficiency as some processors sit idle waiting for the heavily loaded ones to finish.
> + Performing a load balancing check every certain step and see if imbalance has occurred, to simplify, we may dentify the most overloaded and underloaded processors.
> + Then the goal is try to move atoms from the overloaded to underloaded. Update such information and continue simluation.



## Previous work
### pmd.c

This program implements a parallel molecular dynamics simulation, utilizing the **MPI** standard. It provides a detailed record of the molecular dynamics simulation, including:

- Time steps (`stepCount * DeltaT`),
- System temperature (`temperature`),
- Per-atom kinetic energy, potential energy, and total energy (`kinEnergy`, `potEnergy`, `totEnergy`).

Additionally, the program reports the total runtime and communication time (`CPU & COMT`) for performance analysis.

Since this program serves as the foundation for parallel molecular dynamics simulations, all future enhancements will build upon its existing structure and functionality.

### Hilbert and Morton curves

Space-filling curves, notably Hilbert and Morton (Z-order) curves, have been explored for their ability to preserve spatial locality when mapping multi-dimensional data to one-dimensional memory spaces. Previous studies have demonstrated that:

- Hilbert Curves: Offer better locality preservation compared to Morton curves, leading to higher cache hit rates in certain scenarios.
- Morton Curves: Simpler to compute and implement, making them attractive for real-time applications despite slightly lower locality preservation.
- Cache Optimization: Techniques like tiling and blocking have been traditionally used to enhance cache performance, but integrating space-filling curves presents a novel approach.

However, existing research often lacks comprehensive evaluations across diverse applications and does not fully exploit the potential synergy between Hilbert and Morton curves for dynamic cache optimization.

## Proposed Solutions
**Optimization Plan for the Provided Code Using Hilbert and Morton Curves**
1. Implementing Spatial Mapping with Hilbert/Morton Curves
- Assign a **Hilbert** or **Morton** index to each cell to map the 3D space into a 1D sequence.
- Before the force calculation (`compute_accel`), reorder particles in memory based on their cell's curve index.
- Traverse cells in the order defined by the Hilbert or Morton curve. Access neighboring cells in this order for force calculations to leverage the reordered memory layout.

2. Profiling and Performance Metrics
- Compare the execution time, communication cost, and scalability (e.g., strong and weak scaling) of the optimized code with the original version.

**Implement load balancing for distributing computational work more evenly**

Proposed method attempts to maintain an even distribution of atoms across processors by adjusting the number of processors in each dimension. The goal is to reduce load imbalance when the number of atoms per processor deviates significantly from the ideal value. Load balancing occurs at regular intervals.

1. The function first calculates the number of atoms each processor should ideally have based on the total number of atoms and the number of processors.
2. Then it checks the imbalance of atoms assigned to the current processor. If the difference (imbalance) between the number of atoms assigned to this processor and the ideal number  exceeds a certain threshold, it tries to redistribute the atoms across processors. In case of significant imbalance, the topology is adjusted to ensure that atoms are more evenly distributed.

## Expected Results
1. Higher Cache Hit Rates: Reordering particle data based on Hilbert/Morton curves will improve memory locality, leading to significantly reduced cache misses.

2. Reduced Execution Time: Enhanced memory access patterns will optimize the performance of critical sections like compute_accel, resulting in faster force calculations.

3. Improved Scalability: The balance load function periodically checks for these imbalances and redistributes the atoms (or domain) to ensure that each process is responsible for approximately the same number of atoms. This reduces idle time and improves the efficiency of the overall simulation.

4. Reduced Simulation Time: Load balancing can adapt dynamically and can optimize resource utilization, leading to a reduced running time.

## Experiment
### Optimization Plan Using Hilbert and Morton Curves

In our code, the variable cache_accesses keeps track of the total number of cache accesses, while cache_hits tracks how many times data is successfully retrieved from the cache. The cache hit rate can then be calculated by dividing cache_hits by cache_accesses and multiplying by 100 to get the percentage:

$$
\text{Cache Hit Rate} = \frac{\text{Cache Hits}}{\text{Cache Accesses}} \times 100\%
$$

1. Applying Hilbert Curves

| Metric             | Before Using Hilbert Curve | After Using Hilbert Curve |
|--------------------|---------------------------|----------------------------|
| **CPU**            | 5.150883                  | 1.547965                   |
| **COMT**           | 0.02530489                | 0.05067728                 |
| **CATCH HIT RATE** | 14.54%                    | 46.91%                     |

<img src="https://github.com/sc18zs/final_project_cs596/blob/main/Hilbert/Result_comparision.png" alt="Comparison of CPU, COMT, and Cache Hit Rate" width="600"/>

- The Hilbert curve significantly increase the cache hit rate from **14.54%** to the **46.91%**. This substantial improvement indicates that data reordering effectively enhances spatial locality, allowing the simulation to make better use of the cache hierarchy. Higher cache hit rates reduce the frequency of expensive memory accesses, thereby decreasing latency and improving overall performance.

- The Hilbert curve optimization resulted in a **69.95%** reduction in CPU time, demonstrating a significant enhancement in computational efficiency. The slight increase in communication time is minimal and outweighed by the gains in computation speed. This reduction suggests that the optimized memory access patterns allow the processor to execute computations more swiftly by minimizing cache misses and associated delays.

2. Applying Morton Curves

| Metric             | Before Using Morton Curve | After Using Morton Curve |
|--------------------|---------------------------|--------------------------|
| **CPU**            | 0.3005791                 | 0.1117183                |
| **COMT**           | 0.01215395                | 0.01254960               |
| **CATCH HIT RATE** | 15.15%                    | 47.87%                   |

<img src="https://github.com/sc18zs/final_project_cs596/blob/main/Morton/Comparison.png" alt="Comparison of CPU, COMT, and Cache Hit Rate" width="600"/>

After using the Morton Curve, the program's performance showed significant improvement. Specifically, the CPU time decreased by approximately **62.8%**, indicating a substantial increase in computational efficiency and a noticeable speedup of the program. Although COMT time increased by **3.25%**, this change is minimal compared to the reduction in CPU time and can be considered negligible. More importantly, the CATCH HIT RATE increased from **15.15%** to **47.87%**, significantly improving the cache hit rate and reducing memory access failures, further enhancing the program's execution efficiency. Overall, the use of the Morton Curve optimization led to notable improvements in both memory access efficiency and computational performance.

### Applying Load Balancing

**Overview**

This part of report analyzes the implementation of dynamic load balancing functionality in MD.

**Key Changes Overview:**

- Added dynamic load balancing
- Modified process initialization
- Added new configuration parameters
- Modified single_step() structure

**Baseline**

During the simulation, it was observed that significant load imbalances rarely occurred naturally. Please refer to the no imbalance output file under load balancing folder. There are hardly any imbalance at current level.

      === Load Balance Check at Step 50 ===
      Initial atom distribution:
      Process 0: 18473 atoms
      Process 1: 18402 atoms
      Process 2: 18431 atoms
      Process 3: 18470 atoms
      Process 4: 18447 atoms
      Process 5: 18384 atoms
      Process 6: 18386 atoms
      Process 7: 18463 atoms
      Load imbalance metrics:
      Max atoms: 18473, Min atoms: 18384
      Current imbalance ratio: 1.004841
      No load balancing needed at this step.
   

This can be attributed to several key physical and computational factors:
   
   - The molecular dynamics simulation employs a 1×1×2 spatial domain decomposition, where the simulation box is divided only along the z-direction between two processes. The initial configuration uses a face-centered cubic (FCC) lattice structure, which inherently provides a highly uniform distribution of atoms across the domain. This initial setup ensures that each process begins with approximately equal computational load.
   
   
   - During the simulation, it was observed that significant load imbalances rarely occurred naturally. The initial configuration uses a structure, which inherently provides a highly uniform distribution of atoms across the domain. This initial setup ensures that each process begins with approximately equal computational load.
   
   1. Thermodynamic Equilibrium:
      - The system maintains a uniform density due to thermodynamic equilibrium
      - Random thermal motions of atoms tend to preserve the overall spatial distribution
      - Forces between atoms naturally resist clustering or density fluctuations
   
   2. Periodic Boundary Conditions:
      - The implementation of periodic boundary conditions prevents systematic drift
      - When atoms exit one boundary, they enter from the opposite boundary
      - This mechanism maintains a consistent number of atoms in each domain
   
   3. Inter-Process Atom Movement:
      - Atoms moving across process boundaries are typically balanced by reciprocal movements
      - The Lennard-Jones potential maintains relatively uniform density distributions
      - Local density fluctuations tend to be temporary and self-correcting
   

**Load Balance Testing Requirements**

Due to the inherent stability of the atom distribution, artificial imbalance scenarios had to be created to effectively test the load balancing mechanism. The natural dynamics of the system consistently maintained a near-optimal distribution of computational load across processes, with typical imbalance ratios staying well below the triggering threshold of 1.001.


**Mannully Creating a imbalance senario.**

1. Direct Manipulation of Initial Distribution:
    if (sid == 1) {
        // Process 1 keeps only 40% of its atoms
        n = (int)(n * 0.4);
        printf("Process 1: Reduced to %d atoms\n", n);
    } else if (sid == 0) {
        // Process 0 keeps its atoms plus copies some
        for (j = 0; j < original_n * 0.6; j++) {  // Add 60% more
            for (a = 0; a < 3; a++) {
                r[n][a] = r[j][a];  // Copy positions from existing atoms
            }
            n++;
        }
        printf("Process 0: Increased to %d atoms\n", n);
    }

2. System Parameters:
      Reduced the overall system size to make imbalances more noticeable
      Modified the input file (pmdotoc.in) to use smaller InitUcell values
      Adjusted load_imbalance_threshold to ensure the load balancing mechanism would trigger

**Modified Single Step**
We modified single step to handle load balancing in the correct order and without redundant operations. 

**Load Balancing Function**

Key Features of the Load Balancing:
    
    Periodic checks: Runs every load_balance_interval steps
    Threshold-based: Only triggers if imbalance exceeds threshold
    Boundary-aware: Transfers atoms near process boundaries
    Limited transfers: Maximum 10 atoms per step to avoid instability is used in the output in the folder. Large values lead to instability and crash.
    Neighbor-only: Transfers only between adjacent processes
    Preserves physics: Maintains periodic boundary conditions during transfers
     
         === Load Balance Check at Step 100 ===
         loading balancing triggered
         Initial atom distribution:
         Process 0: 11985 atoms
         Process 1: 10705 atoms
         Process 2: 27931 atoms
         Process 3: 26793 atoms
         load threshold is:1.050000
         Load imbalance metrics:
         Max atoms: 27931, Min atoms: 10705
         Current imbalance ratio: 2.609155
         
         Final atom distribution:
         Process 0: 11995 atoms
         Process 1: 10705 atoms
         Process 2: 27921 atoms
         Process 3: 26793 atoms


**Error - numerical overflow**

Extremely large forces happen and indicate a serious issue with particle positions - likely atoms are far too close to each other, causing numerical overflow in the Lennard-Jones force calculation. We could add a minimum distance check and position validation as will be shown in the code. However, this is very computational expensive and not able to get results for serval results. For each atom, we must calculate forces with all nearby atoms. This involves checking atoms in 27 neighboring cells (3×3×3 region)
Complexity is roughly O(N²) in the worst case for N atoms in neighboring cells.

slurmstepd: error: *** JOB 28336171 ON a01-07 CANCELLED AT 2024-12-10T13:05:07 DUE TO TIME LIMIT ***


Square root calculations
Division operations
Multiple floating-point multiplications
Key Changes:


    In compute_accel():
        
        Added minimum separation check (rMin2)
        Added force capping (maxForce)
        Improved cell index calculation with bounds checking
        Added periodic boundary enforcement
    

    In half_kick():

        Added velocity capping (vmax)
        Added warning messages for capped velocities


    In atom_move():

        Added position and velocity validation
        Added checks for NaN and infinity
        Improved periodic boundary enforcement


    In eval_props():

        Added energy validation
        Added error handling for unreasonable energy values
        



## Acknowledgments
This project is part of the CSCI596 course and focuses on enhancing the efficiency of molecular dynamics simulations through innovative computational methods.
