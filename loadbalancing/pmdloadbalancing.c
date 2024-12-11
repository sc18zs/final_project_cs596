/*----------------------------------------------------------------------
Program pmd.c performs parallel molecular-dynamics for Lennard-Jones 
systems using the Message Passing Interface (MPI) standard.
----------------------------------------------------------------------*/
#include "pmdotoc.h"
#include <stdlib.h>
#define min(a,b) ((a) < (b) ? (a) : (b))
/*--------------------------------------------------------------------*/
int main(int argc, char **argv) {
/*--------------------------------------------------------------------*/
  double cpu1;
if (sid == 0) printf("DEBUG: Starting program\n");
  MPI_Init(&argc,&argv); /* Initialize the MPI environment */
    if (sid == 0) printf("DEBUG: After MPI_Init\n");
  MPI_Comm_rank(MPI_COMM_WORLD, &sid);  /* My processor ID */
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (sid == 0) printf("DEBUG: nproc = %d\n", nproc);
    // Set decomposition for 4 processes
    vproc[0] = 1;
    vproc[1] = 1;
    vproc[2] = 4;  // 1×2×2 = 4 processors
    // Verify decomposition
    if (vproc[0] * vproc[1] * vproc[2] != nproc) {
        if (sid == 0) {
            printf("Error: Process grid %d×%d×%d doesn't match nproc=%d\n",
                   vproc[0], vproc[1], vproc[2], nproc);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    
  /* Vector index of this processor */
  vid[0] = sid/(vproc[1]*vproc[2]);
  vid[1] = (sid/vproc[2])%vproc[1];
  vid[2] = sid%vproc[2];


    if (sid == 0) printf("DEBUG: After MPI Init\n");
      
      init_params();
      if (sid == 0) printf("DEBUG: After init_params\n");
      
      set_topology();
      if (sid == 0) printf("DEBUG: After set_topology\n");
      
      init_conf();
      if (sid == 0) printf("DEBUG: After init_conf\n");
      
      if (sid == 0) printf("DEBUG: About to start atom_copy\n");
      atom_copy();
      if (sid == 0) printf("DEBUG: After atom_copy\n");
      
      if (sid == 0) printf("DEBUG: About to start compute_accel\n");
      compute_accel();
      if (sid == 0) printf("DEBUG: After compute_accel\n");

      cpu1 = MPI_Wtime();
      if (sid == 0) printf("DEBUG: Starting main loop\n");
  for (stepCount=1; stepCount<=StepLimit; stepCount++) {
    single_step(); 
    if (stepCount%StepAvg == 0) eval_props();
  }
  cpu = MPI_Wtime() - cpu1;
  if (sid == 0) printf("CPU & COMT = %le %le\n",cpu,comt);

  MPI_Finalize(); /* Clean up the MPI environment */
  return 0;
}

/*--------------------------------------------------------------------*/
void init_params() {
/*----------------------------------------------------------------------
Initializes parameters.
----------------------------------------------------------------------*/
  int a;
  double rr,ri2,ri6,r1;
  FILE *fp;

  /* Read control parameters */
  fp = fopen("pmdotoc.in","r");
  fscanf(fp,"%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
  fscanf(fp,"%le",&Density);
  fscanf(fp,"%le",&InitTemp);
  fscanf(fp,"%le",&DeltaT);
  fscanf(fp,"%d",&StepLimit);
  fscanf(fp,"%d",&StepAvg);
  fscanf(fp,"%d",&load_balance_interval);
  fscanf(fp,"%lf",&load_imbalance_threshold);
  fclose(fp);

  /* Compute basic parameters */
  DeltaTH = 0.5*DeltaT;
  for (a=0; a<3; a++) al[a] = InitUcell[a]/pow(Density/4.0,1.0/3.0);

    printf("Process %d: Initial values:\n", sid);
    printf("Process %d: InitUcell = [%d, %d, %d]\n",
           sid, InitUcell[0], InitUcell[1], InitUcell[2]);
    printf("Process %d: Density = %f\n", sid, Density);
    double pow_term = pow(Density/4.0, 1.0/3.0);
    printf("Process %d: pow(Density/4.0, 1.0/3.0) = %f\n", sid, pow_term);

  if (sid == 0) printf("al = %e %e %e\n",al[0],al[1],al[2]);

  /* Compute the # of cells for linked cell lists */
  for (a=0; a<3; a++) {
    lc[a] = al[a]/RCUT; 
    rc[a] = al[a]/lc[a];
  }
  if (sid == 0) {
    printf("lc = %d %d %d\n",lc[0],lc[1],lc[2]);
    printf("rc = %e %e %e\n",rc[0],rc[1],rc[2]);
  }

  /* Constants for potential truncation */
  rr = RCUT*RCUT; ri2 = 1.0/rr; ri6 = ri2*ri2*ri2; r1=sqrt(rr);
  Uc = 4.0*ri6*(ri6 - 1.0);
  Duc = -48.0*ri6*(ri6 - 0.5)/r1;
}

/*--------------------------------------------------------------------*/
void set_topology() {
/*----------------------------------------------------------------------
Defines a logical network topology.  Prepares a neighbor-node ID table, 
nn, & a shift-vector table, sv, for internode message passing.  Also 
prepares the node parity table, myparity.
----------------------------------------------------------------------*/
  /* Integer vectors to specify the six neighbor nodes */
  int iv[6][3] = {
    {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}
  };
  int ku,a,k1[3];

  /* Set up neighbor tables, nn & sv */
  for (ku=0; ku<6; ku++) {
    /* Vector index of neighbor ku */
    for (a=0; a<3; a++)
      k1[a] = (vid[a]+iv[ku][a]+vproc[a])%vproc[a];
    /* Scalar neighbor ID, nn */
    nn[ku] = k1[0]*vproc[1]*vproc[2]+k1[1]*vproc[2]+k1[2];
    /* Shift vector, sv */
    for (a=0; a<3; a++) sv[ku][a] = al[a]*iv[ku][a];
  }

  /* Set up the node parity table, myparity */
  for (a=0; a<3; a++) {
    if (vproc[a] == 1) 
      myparity[a] = 2;
    else if (vid[a]%2 == 0)
      myparity[a] = 0;
    else
      myparity[a] = 1;
  }
}

/*--------------------------------------------------------------------*/
void init_conf() {
/*----------------------------------------------------------------------
r are initialized to face-centered cubic (fcc) lattice positions.  
rv are initialized with a random velocity corresponding to Temperature.  
----------------------------------------------------------------------*/
  double c[3],gap[3],e[3],vSum[3],gvSum[3],vMag;
  int j,a,nX,nY,nZ;
  double seed;
  /* FCC atoms in the original unit cell */
  double origAtom[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                           {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}}; 

  /* Set up a face-centered cubic (fcc) lattice */
  for (a=0; a<3; a++) gap[a] = al[a]/InitUcell[a];
  n = 0;
  for (nZ=0; nZ<InitUcell[2]; nZ++) {
    c[2] = nZ*gap[2];
    for (nY=0; nY<InitUcell[1]; nY++) {
      c[1] = nY*gap[1];
      for (nX=0; nX<InitUcell[0]; nX++) {
        c[0] = nX*gap[0];
        for (j=0; j<4; j++) {
          for (a=0; a<3; a++)
            r[n][a] = c[a] + gap[a]*origAtom[j][a];
          ++n;
        }
      }
    }
  }
    
    int original_n = n;
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
  /* Total # of atoms summed over processors */
  MPI_Allreduce(&n,&nglob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (sid == 0) printf("nglob = %d\n",nglob);

  /* Generate random velocities */
  seed = 13597.0+sid;
  vMag = sqrt(3*InitTemp);
  for(a=0; a<3; a++) vSum[a] = 0.0;
  for(j=0; j<n; j++) {
    RandVec3(e,&seed);
    for (a=0; a<3; a++) {
      rv[j][a] = vMag*e[a];
      vSum[a] = vSum[a] + rv[j][a];
    }
  }
  MPI_Allreduce(vSum,gvSum,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  /* Make the total momentum zero */
  for (a=0; a<3; a++) gvSum[a] /= nglob;
  for (j=0; j<n; j++)
    for(a=0; a<3; a++) rv[j][a] -= gvSum[a];
}

/*--------------------------------------------------------------------*/
//void single_step() {
///*----------------------------------------------------------------------
//r & rv are propagated by DeltaT using the velocity-Verlet scheme.
//----------------------------------------------------------------------*/
//  int i,a;
//
//  half_kick(); /* First half kick to obtain v(t+Dt/2) */
//  for (i=0; i<n; i++) /* Update atomic coordinates to r(t+Dt) */
//    for (a=0; a<3; a++) r[i][a] = r[i][a] + DeltaT*rv[i][a];
//  atom_move();
//  atom_copy();
//  compute_accel(); /* Computes new accelerations, a(t+Dt) */
//  half_kick(); /* Second half kick to obtain v(t+Dt) */
//}

// In dynamic_load_balance(), remove this line after load balancing
// atom_copy();  // Remove this

// Let the normal single_step() handle atom_copy() instead
void single_step() {
    int i,a;
    
    // Print debug at start
    printf("Process %d: Starting step %d, initial atoms n=%d\n", sid, stepCount, n);
    
    half_kick();
    
    // Update positions
    for (i=0; i<n; i++)
        for (a=0; a<3; a++)
            r[i][a] = r[i][a] + DeltaT*rv[i][a];
            
    printf("Process %d: After position update - n=%d, first atom=(%.3f,%.3f,%.3f)\n",
           sid, n, r[0][0], r[0][1], r[0][2]);
    
    atom_move();
    printf("Process %d: After atom_move - n=%d, first atom=(%.3f,%.3f,%.3f)\n",
           sid, n, r[0][0], r[0][1], r[0][2]);
    
    atom_copy();
    printf("Process %d: After atom_copy - n=%d, nb=%d, first atom=(%.3f,%.3f,%.3f)\n",
           sid, n, nb, r[0][0], r[0][1], r[0][2]);
    
    // Do load balancing after boundaries are handled
    if (stepCount % load_balance_interval == 0) {
        dynamic_load_balance(stepCount);
        printf("Process %d: After load balance - n=%d, first atom=(%.3f,%.3f,%.3f)\n",
               sid, n, r[0][0], r[0][1], r[0][2]);
    }
    
    compute_accel();
    printf("Process %d: After compute_accel - n=%d, first atom force=(%.3f,%.3f,%.3f)\n",
           sid, n, ra[0][0], ra[0][1], ra[0][2]);
           
    half_kick();
    printf("Process %d: End of step - n=%d, first atom vel=(%.3f,%.3f,%.3f)\n",
           sid, n, rv[0][0], rv[0][1], rv[0][2]);
}
/*--------------------------------------------------------------------*/
void half_kick() {
/*----------------------------------------------------------------------
Accelerates atomic velocities, rv, by half the time step.
----------------------------------------------------------------------*/
  int i,a;
  for (i=0; i<n; i++)
    for (a=0; a<3; a++) rv[i][a] = rv[i][a]+DeltaTH*ra[i][a];
}
/*--------------------------------------------------------------------*/
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












void atom_copy() {
/*----------------------------------------------------------------------
Exchanges boundary-atom coordinates among neighbor nodes:  Makes
boundary-atom list, LSB, then sends & receives boundary atoms.
----------------------------------------------------------------------*/
  int kd,kdd,i,ku,inode,nsd,nrc,a;
  int nbnew = 0; /* # of "received" boundary atoms */
  double com1;

/* Main loop over x, y & z directions starts--------------------------*/

  for (kd=0; kd<3; kd++) {

    /* Make a boundary-atom list, LSB---------------------------------*/

    /* Reset the # of to-be-copied atoms for lower&higher directions */
    for (kdd=0; kdd<2; kdd++) lsb[2*kd+kdd][0] = 0;

    /* Scan all the residents & copies to identify boundary atoms */
    for (i=0; i<n+nbnew; i++) {
      for (kdd=0; kdd<2; kdd++) {
        ku = 2*kd+kdd; /* Neighbor ID */
        /* Add an atom to the boundary-atom list, LSB, for neighbor ku
           according to bit-condition function, bbd */
        if (bbd(r[i],ku)) lsb[ku][++(lsb[ku][0])] = i;
      }
    }

    /* Message passing------------------------------------------------*/

    com1=MPI_Wtime(); /* To calculate the communication time */

    /* Loop over the lower & higher directions */
    for (kdd=0; kdd<2; kdd++) {

      inode = nn[ku=2*kd+kdd]; /* Neighbor node ID */

      /* Send & receive the # of boundary atoms-----------------------*/

      nsd = lsb[ku][0]; /* # of atoms to be sent */

      /* Even node: send & recv */
      if (myparity[kd] == 0) {
        MPI_Send(&nsd,1,MPI_INT,inode,10,MPI_COMM_WORLD);
        MPI_Recv(&nrc,1,MPI_INT,MPI_ANY_SOURCE,10,
                 MPI_COMM_WORLD,&status);
      }
      /* Odd node: recv & send */
      else if (myparity[kd] == 1) {
        MPI_Recv(&nrc,1,MPI_INT,MPI_ANY_SOURCE,10,
                 MPI_COMM_WORLD,&status);
        MPI_Send(&nsd,1,MPI_INT,inode,10,MPI_COMM_WORLD);
      }
      /* Single layer: Exchange information with myself */
      else
        nrc = nsd;
      /* Now nrc is the # of atoms to be received */

      /* Send & receive information on boundary atoms-----------------*/

      /* Message buffering */
      for (i=1; i<=nsd; i++)
        for (a=0; a<3; a++) /* Shift the coordinate origin */
          dbuf[3*(i-1)+a] = r[lsb[ku][i]][a]-sv[ku][a];

      /* Even node: send & recv */
      if (myparity[kd] == 0) {
        MPI_Send(dbuf,3*nsd,MPI_DOUBLE,inode,20,MPI_COMM_WORLD);
        MPI_Recv(dbufr,3*nrc,MPI_DOUBLE,MPI_ANY_SOURCE,20,
                 MPI_COMM_WORLD,&status);
      }
      /* Odd node: recv & send */
      else if (myparity[kd] == 1) {
        MPI_Recv(dbufr,3*nrc,MPI_DOUBLE,MPI_ANY_SOURCE,20,
                 MPI_COMM_WORLD,&status);
        MPI_Send(dbuf,3*nsd,MPI_DOUBLE,inode,20,MPI_COMM_WORLD);
      }
      /* Single layer: Exchange information with myself */
      else
        for (i=0; i<3*nrc; i++) dbufr[i] = dbuf[i];

      /* Message storing */
      for (i=0; i<nrc; i++)
        for (a=0; a<3; a++) r[n+nbnew+i][a] = dbufr[3*i+a];

      /* Increment the # of received boundary atoms */
      nbnew = nbnew+nrc;

      /* Internode synchronization */
      MPI_Barrier(MPI_COMM_WORLD);

    } /* Endfor lower & higher directions, kdd */

    comt += MPI_Wtime()-com1; /* Update communication time, COMT */

  } /* Endfor x, y & z directions, kd */

  /* Main loop over x, y & z directions ends--------------------------*/

  /* Update the # of received boundary atoms */
  nb = nbnew;
}



void compute_accel() {
    int i,j,a,lc2[3],lcyz2,lcxyz2,mc[3],c,mc1[3],c1;
    int bintra;
    double dr[3],rr,ri2,ri6,r1,rrCut,fcVal,f,vVal,lpe;

    // Debug at start
    printf("Process %d: compute_accel start, n=%d, nb=%d\n", sid, n, nb);

    /* Reset the potential & forces */
    lpe = 0.0;
    for (i=0; i<n; i++) for (a=0; a<3; a++) ra[i][a] = 0.0;

    /* Make a linked-cell list */
    for (a=0; a<3; a++) lc2[a] = lc[a]+2;
    lcyz2 = lc2[1]*lc2[2];
    lcxyz2 = lc2[0]*lcyz2;

    printf("Process %d: Cell dimensions: lc2=(%d,%d,%d), lcyz2=%d, lcxyz2=%d\n",
           sid, lc2[0], lc2[1], lc2[2], lcyz2, lcxyz2);

    /* Reset the headers */
    for (c=0; c<lcxyz2; c++) head[c] = EMPTY;

    /* Construct headers & linked lists */
    printf("Process %d: Building cell lists for %d total atoms\n", sid, n+nb);
    for (i=0; i<n+nb; i++) {
        // Debug every 1000th atom and first/last few
        if (i < 5 || i > n+nb-5 || i % 1000 == 0) {
            printf("Process %d: Atom %d position: (%.6f,%.6f,%.6f)\n",
                   sid, i, r[i][0], r[i][1], r[i][2]);
        }

        for (a=0; a<3; a++) {
            mc[a] = (r[i][a]+rc[a])/rc[a];
            // Validate cell indices
            if (mc[a] < 0 || mc[a] >= lc2[a]) {
                printf("Process %d: ERROR - Invalid cell index for atom %d: mc=(%d,%d,%d), pos=(%.6f,%.6f,%.6f), rc=(%.6f,%.6f,%.6f)\n",
                       sid, i, mc[0], mc[1], mc[2], r[i][0], r[i][1], r[i][2], rc[0], rc[1], rc[2]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        c = mc[0]*lcyz2+mc[1]*lc2[2]+mc[2];
        
        // Debug cell assignment for specific atoms
        if (i < 5 || i > n+nb-5 || i % 1000 == 0) {
            printf("Process %d: Atom %d assigned to cell %d (mc=%d,%d,%d)\n",
                   sid, i, c, mc[0], mc[1], mc[2]);
        }

        // Validate cell index
        if (c < 0 || c >= lcxyz2) {
            printf("Process %d: ERROR - Invalid cell number %d for atom %d\n", sid, c, i);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        lscl[i] = head[c];
        head[c] = i;
    }

    printf("Process %d: Starting force calculation\n", sid);

    /* Calculate pair interaction */
    rrCut = RCUT*RCUT;

    /* Scan inner cells */
    for (mc[0]=1; mc[0]<=lc[0]; (mc[0])++)
    for (mc[1]=1; mc[1]<=lc[1]; (mc[1])++)
    for (mc[2]=1; mc[2]<=lc[2]; (mc[2])++) {
        /* Calculate a scalar cell index */
        c = mc[0]*lcyz2+mc[1]*lc2[2]+mc[2];
        
        /* Skip this cell if empty */
        if (head[c] == EMPTY) continue;

        /* Scan the neighbor cells (including itself) of cell c */
        for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
        for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++)
        for (mc1[2]=mc[2]-1; mc1[2]<=mc[2]+1; (mc1[2])++) {
            /* Calculate the scalar cell index of the neighbor cell */
            c1 = mc1[0]*lcyz2+mc1[1]*lc2[2]+mc1[2];
            /* Skip this neighbor cell if empty */
            if (head[c1] == EMPTY) continue;

            /* Scan atom i in cell c */
            i = head[c];
            while (i != EMPTY) {
                /* Scan atom j in cell c1 */
                j = head[c1];
                while (j != EMPTY) {
                    /* No calculation with itself */
                    if (j != i) {
                        /* Logical flag: intra(true)- or inter(false)-pair atom */
                        bintra = (j < n);

                        /* Pair vector dr = r[i] - r[j] */
                        for (rr=0.0, a=0; a<3; a++) {
                            dr[a] = r[i][a]-r[j][a];
                            rr += dr[a]*dr[a];
                        }

                        /* Calculate forces if distance between atoms is less than cut-off */
                        if (i<j && rr<rrCut) {
                            ri2 = 1.0/rr;
                            ri6 = ri2*ri2*ri2;
                            r1 = sqrt(rr);
                            fcVal = 48.0*ri2*ri6*(ri6-0.5) + Duc/r1;
                            vVal = 4.0*ri6*(ri6-1.0) - Uc - Duc*(r1-RCUT);

                            // Debug any suspiciously large forces
                            if (fabs(fcVal) > 1000.0) {
                                printf("Process %d: WARNING - Large force between atoms %d and %d: %.6f\n",
                                       sid, i, j, fcVal);
                                printf("Process %d: Positions: (%.6f,%.6f,%.6f) and (%.6f,%.6f,%.6f)\n",
                                       sid, r[i][0], r[i][1], r[i][2], r[j][0], r[j][1], r[j][2]);
                            }

                            if (bintra) lpe += vVal;
                            else lpe += 0.5*vVal;

                            for (a=0; a<3; a++) {
                                f = fcVal*dr[a];
                                ra[i][a] += f;
                                if (bintra) ra[j][a] -= f;
                            }
                        }
                    }
                    j = lscl[j];
                }
                i = lscl[i];
            }
        }
    }

    /* Global potential energy */
    MPI_Allreduce(&lpe, &potEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    printf("Process %d: compute_accel complete, potEnergy=%.6f\n", sid, potEnergy);
}

/*--------------------------------------------------------------------*/
void eval_props() {
/*----------------------------------------------------------------------
Evaluates physical properties: kinetic, potential & total energies.
----------------------------------------------------------------------*/
  double vv,lke;
  int i,a;

  /* Total kinetic energy */
  for (lke=0.0, i=0; i<n; i++) {
    for (vv=0.0, a=0; a<3; a++) vv += rv[i][a]*rv[i][a];
    lke += vv;
  }
  lke *= 0.5;
  MPI_Allreduce(&lke,&kinEnergy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  /* Energy paer atom */
  kinEnergy /= nglob;
  potEnergy /= nglob;
  totEnergy = kinEnergy + potEnergy;
  temperature = kinEnergy*2.0/3.0;

  /* Print the computed properties */
  if (sid == 0) printf("Computed Properties %9.6f %9.6f %9.6f %9.6f\n",
                stepCount*DeltaT,temperature,potEnergy,totEnergy);
}

/*--------------------------------------------------------------------*/
void atom_move() {
/*----------------------------------------------------------------------
Sends moved-out atoms to neighbor nodes and receives moved-in atoms 
from neighbor nodes.  Called with n, r[0:n-1] & rv[0:n-1], atom_move 
returns a new n' together with r[0:n'-1] & rv[0:n'-1].
----------------------------------------------------------------------*/

/* Local variables------------------------------------------------------

mvque[6][NBMAX]: mvque[ku][0] is the # of to-be-moved atoms to neighbor 
  ku; MVQUE[ku][k>0] is the atom ID, used in r, of the k-th atom to be
  moved.
----------------------------------------------------------------------*/
  int mvque[6][NBMAX];
  int newim = 0; /* # of new immigrants */
  int ku,kd,i,kdd,kul,kuh,inode,ipt,a,nsd,nrc;
  double com1;

  /* Reset the # of to-be-moved atoms, MVQUE[][0] */
  for (ku=0; ku<6; ku++) mvque[ku][0] = 0;

  /* Main loop over x, y & z directions starts------------------------*/

  for (kd=0; kd<3; kd++) {

    /* Make a moved-atom list, mvque----------------------------------*/

    /* Scan all the residents & immigrants to list moved-out atoms */
    for (i=0; i<n+newim; i++) {
      kul = 2*kd  ; /* Neighbor ID */
      kuh = 2*kd+1; 
      /* Register a to-be-copied atom in mvque[kul|kuh][] */      
      if (r[i][0] > MOVED_OUT) { /* Don't scan moved-out atoms */
        /* Move to the lower direction */
        if (bmv(r[i],kul)) mvque[kul][++(mvque[kul][0])] = i;
        /* Move to the higher direction */
        else if (bmv(r[i],kuh)) mvque[kuh][++(mvque[kuh][0])] = i;
      }
    }

    /* Message passing with neighbor nodes----------------------------*/

    com1 = MPI_Wtime();

    /* Loop over the lower & higher directions------------------------*/

    for (kdd=0; kdd<2; kdd++) {

      inode = nn[ku=2*kd+kdd]; /* Neighbor node ID */

      /* Send atom-number information---------------------------------*/  

      nsd = mvque[ku][0]; /* # of atoms to-be-sent */

      /* Even node: send & recv */
      if (myparity[kd] == 0) {
        MPI_Send(&nsd,1,MPI_INT,inode,110,MPI_COMM_WORLD);
        MPI_Recv(&nrc,1,MPI_INT,MPI_ANY_SOURCE,110,
                 MPI_COMM_WORLD,&status);
      }
      /* Odd node: recv & send */
      else if (myparity[kd] == 1) {
        MPI_Recv(&nrc,1,MPI_INT,MPI_ANY_SOURCE,110,
                 MPI_COMM_WORLD,&status);
        MPI_Send(&nsd,1,MPI_INT,inode,110,MPI_COMM_WORLD);
      }
      /* Single layer: Exchange information with myself */
      else
        nrc = nsd;
      /* Now nrc is the # of atoms to be received */

      /* Send & receive information on boundary atoms-----------------*/

      /* Message buffering */
      for (i=1; i<=nsd; i++)
        for (a=0; a<3; a++) {
          /* Shift the coordinate origin */
          dbuf[6*(i-1)  +a] = r [mvque[ku][i]][a]-sv[ku][a]; 
          dbuf[6*(i-1)+3+a] = rv[mvque[ku][i]][a];
          r[mvque[ku][i]][0] = MOVED_OUT; /* Mark the moved-out atom */
        }

      /* Even node: send & recv, if not empty */
      if (myparity[kd] == 0) {
        MPI_Send(dbuf,6*nsd,MPI_DOUBLE,inode,120,MPI_COMM_WORLD);
        MPI_Recv(dbufr,6*nrc,MPI_DOUBLE,MPI_ANY_SOURCE,120,
                 MPI_COMM_WORLD,&status);
      }
      /* Odd node: recv & send, if not empty */
      else if (myparity[kd] == 1) {
        MPI_Recv(dbufr,6*nrc,MPI_DOUBLE,MPI_ANY_SOURCE,120,
                 MPI_COMM_WORLD,&status);
        MPI_Send(dbuf,6*nsd,MPI_DOUBLE,inode,120,MPI_COMM_WORLD);
      }
      /* Single layer: Exchange information with myself */
      else
        for (i=0; i<6*nrc; i++) dbufr[i] = dbuf[i];

      /* Message storing */
      for (i=0; i<nrc; i++)
        for (a=0; a<3; a++) {
          r [n+newim+i][a] = dbufr[6*i  +a]; 
          rv[n+newim+i][a] = dbufr[6*i+3+a];
        }

      /* Increment the # of new immigrants */
      newim = newim+nrc;

      /* Internode synchronization */
      MPI_Barrier(MPI_COMM_WORLD);

    } /* Endfor lower & higher directions, kdd */

    comt=comt+MPI_Wtime()-com1;

  } /* Endfor x, y & z directions, kd */
  
  /* Main loop over x, y & z directions ends--------------------------*/

  /* Compress resident arrays including new immigrants */

  ipt = 0;
  for (i=0; i<n+newim; i++) {
    if (r[i][0] > MOVED_OUT) {
      for (a=0; a<3; a++) {
        r [ipt][a] = r [i][a];
        rv[ipt][a] = rv[i][a];
      }
      ++ipt;
    }
  }

  /* Update the compressed # of resident atoms */
  n = ipt;
}

/*----------------------------------------------------------------------
Bit condition functions:

1. bbd(ri,ku) is .true. if coordinate ri[3] is in the boundary to 
     neighbor ku.
2. bmv(ri,ku) is .true. if an atom with coordinate ri[3] has moved out 
     to neighbor ku.
----------------------------------------------------------------------*/
int bbd(double* ri, int ku) {
  int kd,kdd;
  kd = ku/2; /* x(0)|y(1)|z(2) direction */
  kdd = ku%2; /* Lower(0)|higher(1) direction */
  if (kdd == 0)
    return ri[kd] < RCUT;
  else
    return al[kd]-RCUT < ri[kd];
}
int bmv(double* ri, int ku) {
  int kd,kdd;
  kd = ku/2; /* x(0)|y(1)|z(2) direction */
  kdd = ku%2; /* Lower(0)|higher(1) direction */
  if (kdd == 0)
    return ri[kd] < 0.0;
  else
    return al[kd] < ri[kd];
}
