/*----------------------------------------------------------------------
Program pmd.c performs parallel molecular-dynamics for Lennard-Jones 
systems using the Message Passing Interface (MPI) standard.
----------------------------------------------------------------------*/
#include "test_pmd.h"

/* Global variables for cache hit tracking */
int cache_hits = 0;  // Count of cache hits
int cache_accesses = 0;  // Total cache accesses



/*--------------------------------------------------------------------*/
int main(int argc, char **argv) {
/*--------------------------------------------------------------------*/
  double cpu1;

  MPI_Init(&argc,&argv); /* Initialize the MPI environment */
  MPI_Comm_rank(MPI_COMM_WORLD, &sid);  /* My processor ID */
  /* Vector index of this processor */
  vid[0] = sid/(vproc[1]*vproc[2]);
  vid[1] = (sid/vproc[2])%vproc[1];
  vid[2] = sid%vproc[2];

  init_params();
  set_topology();
  init_conf();
  assign_morton_indices();
  sort_particles_by_morton();
  atom_copy();
  compute_accel(); /* Computes initial accelerations */

  cpu1 = MPI_Wtime();
  for (stepCount=1; stepCount<=StepLimit; stepCount++) {
    single_step(); 
    if (stepCount%StepAvg == 0) eval_props();
  }
  cpu = MPI_Wtime() - cpu1;
  if (sid == 0) printf("CPU & COMT = %le %le\n",cpu,comt);
    
/* Output cache hit rate */
  if (sid == 0) {
    double hit_rate = (cache_accesses > 0) ? (double)cache_hits / cache_accesses : 0.0;
    printf("Cache hit rate: %.2f%%\n", hit_rate * 100.0);
  }

//  free_morton_indices();
  MPI_Finalize(); /* Clean up the MPI environment */
  return 0;
}

/*将三维坐标通过morton curve映射为一维index*/
int compute_morton_index(int x, int y, int z) {
    if (x >= MAX_MORTON_RANGE || y >= MAX_MORTON_RANGE || z >= MAX_MORTON_RANGE) {
        printf("Error: Morton index out of bounds (x=%d, y=%d, z=%d). Max range = %d\n", x, y, z, MAX_MORTON_RANGE);
        exit(EXIT_FAILURE);
    }

    int morton_index = 0;
    for (int i = 0; i < (sizeof(int) * 8); i++) {
        if ((1 << i) >= MAX_MORTON_RANGE) break; // 确保 Morton 索引位数不超过范围
        morton_index |= ((x >> i & 1) << (3 * i)) |
                        ((y >> i & 1) << (3 * i + 1)) |
                        ((z >> i & 1) << (3 * i + 2));
    }
    return morton_index;
}



/* 将morton index用于粒子映射 */
void assign_morton_indices() {
    for (int i = 0; i < n; i++) {
        int x = (int)(r[i][0] / rc[0]);  // 粒子所在的网格 x 坐标
        int y = (int)(r[i][1] / rc[1]);  // 粒子所在的网格 y 坐标
        int z = (int)(r[i][2] / rc[2]);  // 粒子所在的网格 z 坐标

        if (x < 0 || x >= lc[0] || y < 0 || y >= lc[1] || z < 0 || z >= lc[2]) {
            printf("Error: Particle at index %d has out-of-bounds coordinates (x=%d, y=%d, z=%d).\n", i, x, y, z);
            exit(EXIT_FAILURE);
        }

        particle_indices[i] = cell_indices[x][y][z];  // 映射粒子到 Morton Index
    }
}

/* 确保粒子按照 morton index 进行排序 */
void sort_particles_by_morton() {
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (particle_indices[j] > particle_indices[j + 1]) {
                // Swap Morton indices
                int temp_idx = particle_indices[j];
                particle_indices[j] = particle_indices[j + 1];
                particle_indices[j + 1] = temp_idx;

                // Swap particle positions
                double temp_r[3], temp_v[3];
                for (int k = 0; k < 3; k++) {
                    temp_r[k] = r[j][k];
                    r[j][k] = r[j + 1][k];
                    r[j + 1][k] = temp_r[k];

                    temp_v[k] = rv[j][k];
                    rv[j][k] = rv[j + 1][k];
                    rv[j + 1][k] = temp_v[k];
                }
            }
        }
    }
}

void reorganize_particle_data() {
    // Temporary arrays to store reordered particle data
    double temp_r[n][3];
    double temp_v[n][3];

    // Reorganize based on Morton indices
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            temp_r[i][j] = r[particle_indices[i]][j];
            temp_v[i][j] = rv[particle_indices[i]][j];
        }
    }

    // Copy reordered data back to original arrays
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++) {
            r[i][j] = temp_r[i][j];
            rv[i][j] = temp_v[i][j];
        }
    }
}

void free_morton_indices() {
    for (int x = 0; x < MAX_MORTON_RANGE; x++) {
        for (int y = 0; y < MAX_MORTON_RANGE; y++) {
            free(cell_indices[x][y]);
        }
        free(cell_indices[x]);
    }
    free(cell_indices);
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
  fp = fopen("test_pmd.in","r");
  fscanf(fp,"%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
  fscanf(fp,"%le",&Density);
  fscanf(fp,"%le",&InitTemp);
  fscanf(fp,"%le",&DeltaT);
  fscanf(fp,"%d",&StepLimit);
  fscanf(fp,"%d",&StepAvg);
  fclose(fp);

  /* Compute basic parameters */
  DeltaTH = 0.5*DeltaT;
  for (a=0; a<3; a++) al[a] = InitUcell[a]/pow(Density/4.0,1.0/3.0);
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
    
    cell_indices = (int ***)malloc(lc[0] * sizeof(int **));
      if (cell_indices == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for cell_indices.\n");
        exit(EXIT_FAILURE);
      }

      for (int x = 0; x < lc[0]; x++) {
        cell_indices[x] = (int **)malloc(lc[1] * sizeof(int *));
        if (cell_indices[x] == NULL) {
          fprintf(stderr, "Error: Memory allocation failed for cell_indices[%d].\n", x);
          exit(EXIT_FAILURE);
        }
        for (int y = 0; y < lc[1]; y++) {
          cell_indices[x][y] = (int *)malloc(lc[2] * sizeof(int));
          if (cell_indices[x][y] == NULL) {
            fprintf(stderr, "Error: Memory allocation failed for cell_indices[%d][%d].\n", x, y);
            exit(EXIT_FAILURE);
          }
          for (int z = 0; z < lc[2]; z++) {
            cell_indices[x][y][z] = compute_morton_index(x, y, z);
          }
        }
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
    
  // 输出当前进程的 n 值
  printf("Process %d: n = %d\n", sid, n);
    
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

void single_step() {
    half_kick(); /* First half kick to obtain v(t+Dt/2) */
    for (int i = 0; i < n; i++) /* Update atomic coordinates to r(t+Dt) */
        for (int a = 0; a < 3; a++) r[i][a] = r[i][a] + DeltaT * rv[i][a];
    atom_move();
    assign_morton_indices();
    sort_particles_by_morton();
    reorganize_particle_data();
    atom_copy();
    compute_accel(); /* Computes new accelerations, a(t+Dt) */
    half_kick(); /* Second half kick to obtain v(t+Dt) */
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
/*----------------------------------------------------------------------
Given atomic coordinates, r[0:n+nb-1][], for the extended (i.e.,
resident & copied) system, computes the acceleration, ra[0:n-1][], for
the residents. Morton curve optimization is applied to limit neighbor
search to particles with adjacent Morton indices.
----------------------------------------------------------------------*/
  int i, j, a;
  double dr[3], rr, ri2, ri6, r1, rrCut, fcVal, f, vVal, lpe;

  /* Reset the potential & forces */
  lpe = 0.0;
  for (i = 0; i < n; i++) {
    for (a = 0; a < 3; a++) ra[i][a] = 0.0;
  }

  /* Cut-off distance squared */
  rrCut = RCUT * RCUT;

  /* Compute pair interactions using Morton indices */
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
        /* Increment total cache accesses */
        cache_accesses++;
        
      /* Stop checking if Morton indices are too far apart */
      if (particle_indices[j] - particle_indices[i] > MAX_MORTON_RANGE) {
        break;
      }

      /* Pair vector dr = r[i] - r[j] */
      for (rr = 0.0, a = 0; a < 3; a++) {
        dr[a] = r[i][a] - r[j][a];
        rr += dr[a] * dr[a];
      }

      /* Apply cut-off */
      if (rr < rrCut) {
          /* Increment cache hits */
          cache_hits++;
          
        ri2 = 1.0 / rr;
        ri6 = ri2 * ri2 * ri2;
        r1 = sqrt(rr);
        fcVal = 48.0 * ri2 * ri6 * (ri6 - 0.5) + Duc / r1;
        vVal = 4.0 * ri6 * (ri6 - 1.0) - Uc - Duc * (r1 - RCUT);

        lpe += vVal;
        for (a = 0; a < 3; a++) {
          f = fcVal * dr[a];
          ra[i][a] += f;
          ra[j][a] -= f;
        }
      }
    }
  }

  /* Global potential energy */
  MPI_Allreduce(&lpe, &potEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    printf("Cache accesses: %d, Cache hits: %d, Cache hit rate: %.2f%%\n",
           cache_accesses, cache_hits, (double)cache_hits / cache_accesses * 100);

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
  if (sid == 0) printf("%9.6f %9.6f %9.6f %9.6f\n",
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
        
    /* Ensure particles stay within the box by wrapping around */
          for (a = 0; a < 3; a++) {
            while (r[i][a] < 0.0) r[i][a] += al[a];
            while (r[i][a] >= al[a]) r[i][a] -= al[a];
          }
        
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
