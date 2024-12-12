**Overview**

This report analyzes the implementation of dynamic load balancing functionality in a parallel molecular dynamics simulation program using MPI (Message Passing Interface). The original code was enhanced to handle uneven particle distributions across processes.

**Key Changes Overview:**

- Added dynamic load balancing
- Modified process initialization
- Enhanced debugging/monitoring
- Added new configuration parameters
- Modified single_step() structure

**Baseline**
Please refer to the no imbalance output file. There are hardly any imbalance at current level. 

During the simulation, it was observed that significant load imbalances rarely occurred naturally. This can be attributed to several key physical and computational factors:

The molecular dynamics simulation employs a 1×1×2 spatial domain decomposition, where the simulation box is divided only along the z-direction between two processes. The initial configuration uses a face-centered cubic (FCC) lattice structure, which inherently provides a highly uniform distribution of atoms across the domain. This initial setup ensures that each process begins with approximately equal computational load.


During the simulation, it was observed that significant load imbalances rarely occurred naturally. The initial configuration uses a structure, which inherently provides a highly uniform distribution of atoms across the domain. This initial setup ensures that each process begins with approximately equal computational load.This can be attributed to several key physical and computational factors:

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


Load Balance Testing Requirements

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
3. 
Reduced the overall system size to make imbalances more noticeable
Modified the input file (pmdotoc.in) to use smaller InitUcell values
Adjusted load_imbalance_threshold to ensure the load balancing mechanism would trigger

**Modified Single Step**


**Load Balancing Function**

Key Features of the Load Balancing:
    
    Periodic checks: Runs every load_balance_interval steps
    Threshold-based: Only triggers if imbalance exceeds threshold
    Boundary-aware: Transfers atoms near process boundaries
    Limited transfers: Maximum 10 atoms per step to avoid instability is used in the output in the folder. Large values lead to instability and crash.
    Neighbor-only: Transfers only between adjacent processes
    Preserves physics: Maintains periodic boundary conditions during transfers


    


**=== Load Balance Check at Step 100 ===**
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
        



**Code**


Key Components
Load Imbalance Detection

Gathers atom counts from all processes using MPI_Allgather
Calculates imbalance ratio: max_atoms/min_atoms
Triggers rebalancing if ratio exceeds load_imbalance_threshold

Load Balancing Strategy

Uses a spatial-based approach for 1×1×2 decomposition
Transfers atoms between direct neighbors only (process 0 ↔ 1)
Limits transfers to 1000 atoms per step (max_transfer)
Targets atoms near the boundary between processes:

Process 0: transfers atoms with z > (al[2]/2.0 - RCUT)
Process 1: transfers atoms with z < RCUT



void dynamic_load_balance(int step) {
    int total_procs, max_atoms, min_atoms;
    int *local_atom_counts;
    int *send_counts, *recv_counts;
    int *send_displs, *recv_displs;
    double *send_buffer = NULL, *recv_buffer = NULL;
    int i, j, a;

    MPI_Comm_size(MPI_COMM_WORLD, &total_procs);

    // Allocate arrays
    local_atom_counts = (int*)malloc(total_procs * sizeof(int));
    send_counts = (int*)malloc(total_procs * sizeof(int));
    recv_counts = (int*)malloc(total_procs * sizeof(int));
    send_displs = (int*)malloc(total_procs * sizeof(int));
    recv_displs = (int*)malloc(total_procs * sizeof(int));

    // Check allocations
    if (!local_atom_counts || !send_counts || !recv_counts ||
        !send_displs || !recv_displs) {
        printf("Process %d: Memory allocation failed in load balancing\n", sid);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Debug: Print initial distribution
    if (sid == 0) printf("\n=== Load Balance Check at Step %d ===\n", step);
    
    // Gather current atom distribution
    MPI_Allgather(&n, 1, MPI_INT, local_atom_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
    if (sid == 0) {
        printf("Initial atom distribution:\n");
        for (i = 0; i < total_procs; i++) {
            printf("Process %d: %d atoms\n", i, local_atom_counts[i]);
        }
    }

    // Calculate imbalance
    max_atoms = min_atoms = local_atom_counts[0];
    for (i = 1; i < total_procs; i++) {
        if (local_atom_counts[i] > max_atoms) max_atoms = local_atom_counts[i];
        if (local_atom_counts[i] < min_atoms) min_atoms = local_atom_counts[i];
    }

    if (sid == 0) {
        printf("load threshold is:%f\n", load_imbalance_threshold);
        printf("Load imbalance metrics:\n");
        printf("Max atoms: %d, Min atoms: %d\n", max_atoms, min_atoms);
        printf("Current imbalance ratio: %f\n", (double)max_atoms/min_atoms);
    }

    // Check if load balancing is needed
    if ((double)max_atoms / min_atoms > load_imbalance_threshold) {
        if (sid == 0) printf("\nLoad balancing triggered!\n");

        // Calculate target distribution
        int target_atoms = nglob / total_procs;
        int my_excess = local_atom_counts[sid] - target_atoms;
       // int max_transfer = (local_atom_counts[sid] * 5) / 100; // 5% limit
        int max_transfer = 1000;  // Only transfer up to 1000 atoms at a time
        // Initialize transfer arrays
        for (i = 0; i < total_procs; i++) {
            send_counts[i] = 0;
            recv_counts[i] = 0;
        }

        // Determine transfers if we have excess atoms
        // In dynamic_load_balance
        // Determine transfers if we have excess atoms
        if (my_excess > 0) {
            // For 1×1×2 decomposition, only transfer to direct neighbor in z
            int partner = (sid == 0) ? 1 : 0;
            int their_excess = local_atom_counts[partner] - target_atoms;
            
            // Initialize send_counts
            for (i = 0; i < total_procs; i++) {
                send_counts[i] = 0;
            }
            
            if (their_excess < 0) {
                // We have excess, they have deficit
                int transfer = min(my_excess / 2, abs(their_excess));
                transfer = min(transfer, max_transfer);
                
                // Only transfer atoms near the boundary with neighbor
                if (sid == 0) {
                    // Process 0: transfer atoms with largest z coordinates
                    int count = 0;
                    for (i = 0; i < n && count < transfer; i++) {
                        if (r[i][2] > al[2]/2.0 - RCUT) {  // Near upper boundary
                            // Swap this atom to the transfer region
                            if (i < n - count - 1) {
                                for (a = 0; a < 3; a++) {
                                    double tmp_r = r[i][a];
                                    double tmp_rv = rv[i][a];
                                    r[i][a] = r[n-count-1][a];
                                    rv[i][a] = rv[n-count-1][a];
                                    r[n-count-1][a] = tmp_r;
                                    rv[n-count-1][a] = tmp_rv;
                                }
                            }
                            count++;
                        }
                    }
                    send_counts[partner] = count;
                } else {
                    // Process 1: transfer atoms with smallest z coordinates
                    int count = 0;
                    for (i = 0; i < n && count < transfer; i++) {
                        if (r[i][2] < RCUT) {  // Near lower boundary
                            // Swap this atom to the transfer region
                            if (i < n - count - 1) {
                                for (a = 0; a < 3; a++) {
                                    double tmp_r = r[i][a];
                                    double tmp_rv = rv[i][a];
                                    r[i][a] = r[n-count-1][a];
                                    rv[i][a] = rv[n-count-1][a];
                                    r[n-count-1][a] = tmp_r;
                                    rv[n-count-1][a] = tmp_rv;
                                }
                            }
                            count++;
                        }
                    }
                    send_counts[partner] = count;
                }
            }
        }
        
        
        
        
        
        // Exchange transfer counts
        for (i = 0; i < total_procs; i++) {
            if (i == sid) continue;
            int partner_send_count;
            MPI_Sendrecv(&send_counts[i], 1, MPI_INT, i, 0,
                        &partner_send_count, 1, MPI_INT, i, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recv_counts[i] = partner_send_count;
        }

        // Calculate total transfers and displacements
        int total_send = 0;
        int total_recv = 0;
        send_displs[0] = 0;
        recv_displs[0] = 0;

        for (i = 0; i < total_procs; i++) {
            total_send += send_counts[i];
            total_recv += recv_counts[i];
            if (i > 0) {
                send_displs[i] = send_displs[i-1] + send_counts[i-1];
                recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
            }
        }

        // Debug output
        printf("Process %d: Total atoms to send=%d, receive=%d\n",
               sid, total_send, total_recv);

        // Allocate transfer buffers if needed
        if (total_send > 0) {
            send_buffer = (double*)malloc(total_send * 6 * sizeof(double));
            if (!send_buffer) {
                printf("Process %d: Failed to allocate send buffer\n", sid);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        if (total_recv > 0) {
            recv_buffer = (double*)malloc(total_recv * 6 * sizeof(double));
            if (!recv_buffer) {
                if (send_buffer) free(send_buffer);
                printf("Process %d: Failed to allocate receive buffer\n", sid);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Pack send buffer if we're sending atoms
        if (total_send > 0) {
            int idx = 0;
            int atom_start = n - total_send;
            
            printf("Process %d: Packing atoms from index %d to %d\n",
                   sid, atom_start, n-1);
            
            for (i = 0; i < total_send; i++) {
                int src_idx = atom_start + i;
                
                // Apply periodic boundaries before packing
                for (a = 0; a < 3; a++) {
                    while (r[src_idx][a] < 0) r[src_idx][a] += al[a];
                    while (r[src_idx][a] >= al[a]) r[src_idx][a] -= al[a];
                }
                
                // Print values before packing
                printf("Process %d: Pre-pack atom %d: pos=(%.3f,%.3f,%.3f) vel=(%.3f,%.3f,%.3f)\n",
                       sid, src_idx,
                       r[src_idx][0], r[src_idx][1], r[src_idx][2],
                       rv[src_idx][0], rv[src_idx][1], rv[src_idx][2]);
                       
                // Pack positions
                send_buffer[idx++] = r[src_idx][0];
                send_buffer[idx++] = r[src_idx][1];
                send_buffer[idx++] = r[src_idx][2];
                // Pack velocities
                send_buffer[idx++] = rv[src_idx][0];
                send_buffer[idx++] = rv[src_idx][1];
                send_buffer[idx++] = rv[src_idx][2];
            }
            
            // Print first packed atom
            printf("Process %d: First packed values: (%.3f,%.3f,%.3f)(%.3f,%.3f,%.3f)\n",
                   sid, send_buffer[0], send_buffer[1], send_buffer[2],
                   send_buffer[3], send_buffer[4], send_buffer[5]);
                   
            n -= total_send;
        }

        // Handle MPI communication
        for (i = 0; i < total_procs; i++) {
            send_counts[i] *= 6;
            recv_counts[i] *= 6;
            send_displs[i] *= 6;
            recv_displs[i] *= 6;
            printf("Process %d: Rank %d - send=%d recv=%d sdispl=%d rdispl=%d\n",
                   sid, i, send_counts[i], recv_counts[i], send_displs[i], recv_displs[i]);
        }

        int err = MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_DOUBLE,
                                recv_buffer, recv_counts, recv_displs, MPI_DOUBLE,
                                MPI_COMM_WORLD);

        for (i = 0; i < total_procs; i++) {
            send_counts[i] /= 6;
            recv_counts[i] /= 6;
            send_displs[i] /= 6;
            recv_displs[i] /= 6;
        }

        // Unpack received atoms
        if (total_recv > 0) {
            int idx = 0;
            int old_n = n;
            
            printf("Process %d: Unpacking %d atoms starting at index %d\n",
                   sid, total_recv, n);
            
            // Print first received atom before unpacking
            if (total_recv > 0) {
                printf("Process %d: First received values: (%.3f,%.3f,%.3f)(%.3f,%.3f,%.3f)\n",
                       sid, recv_buffer[0], recv_buffer[1], recv_buffer[2],
                       recv_buffer[3], recv_buffer[4], recv_buffer[5]);
            }
            
            n += total_recv;
            
            for (i = 0; i < total_recv; i++) {
                int dest_idx = old_n + i;
                // Unpack positions
                r[dest_idx][0] = recv_buffer[idx++];
                r[dest_idx][1] = recv_buffer[idx++];
                r[dest_idx][2] = recv_buffer[idx++];
                // Unpack velocities
                rv[dest_idx][0] = recv_buffer[idx++];
                rv[dest_idx][1] = recv_buffer[idx++];
                rv[dest_idx][2] = recv_buffer[idx++];
                
                // Print after unpacking
                printf("Process %d: Post-unpack atom %d: pos=(%.3f,%.3f,%.3f) vel=(%.3f,%.3f,%.3f)\n",
                       sid, dest_idx,
                       r[dest_idx][0], r[dest_idx][1], r[dest_idx][2],
                       rv[dest_idx][0], rv[dest_idx][1], rv[dest_idx][2]);
            }
        }

        // Clean up transfer buffers
        if (send_buffer) free(send_buffer);
        if (recv_buffer) free(recv_buffer);

        // Update neighbor lists and boundary atoms
        atom_copy();

        // Verify new distribution
        MPI_Allgather(&n, 1, MPI_INT, local_atom_counts, 1, MPI_INT, MPI_COMM_WORLD);
        if (sid == 0) {
            printf("\nFinal atom distribution:\n");
            for (i = 0; i < total_procs; i++) {
                printf("Process %d: %d atoms\n", i, local_atom_counts[i]);
            }
        }
    } else {
        if (sid == 0) printf("No load balancing needed at this step.\n");
    }

    // Clean up arrays
    free(local_atom_counts);
    free(send_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);
}
