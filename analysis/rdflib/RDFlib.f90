
!!! This routine calculates all 3 of monomer-monomer, water-water
!!! water-monomer rdfs. It explicitly assumes the following typesettings:
!!! hydrogen - 1, oxygen - 2, monomer - 3. The rdf s are ordered as:
!!! monomer-monomer: 1, water-water : 2, water-monomer: 3

SUBROUTINE find_RDF(g, Nbins, Bin_centers, Bin_delta, &
AtomTypes, Pos, NAtom, BoxL, atomtypedefs)
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: NAtom, Nbins
	REAL(8), INTENT(IN) :: Bin_delta
	REAL(8), INTENT(IN), DIMENSION(0:Nbins-1) :: Bin_centers
	REAL(8), INTENT(IN) :: BoxL
	REAL(8), INTENT(IN), DIMENSION(0:NAtom-1,0:2) :: Pos
	INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
	INTEGER, INTENT(IN), DIMENSION(0:2) :: atomtypedefs
	
	REAL(8), INTENT(INOUT), DIMENSION(0:Nbins-1, 0:2) :: g
	!f2py intent(in,out) :: g
	
	INTEGER :: i,j, Bin_assignment, rdf_type = -1
	REAL(8) :: r, rsq, Bin_min, Bin_max, Bin_minsq, Bin_maxsq, invBoxL, invBin_delta
	REAL(8), DIMENSION(0:2) :: Posi, rij
	INTEGER :: monomer_atomtype, water_atomtype, hydrogen_atomtype
	
	! Precompute bin extrema for fast comparison
	Bin_min  = MINVAL(Bin_centers) - Bin_delta/2.0
	Bin_max = MAXVAL(Bin_centers) + Bin_delta/2.0
	Bin_minsq = Bin_min * Bin_min
	Bin_maxsq = Bin_max * Bin_max
	invBoxL = 1.d0/BoxL
	invBin_delta = 1.d0/Bin_delta
	
	! define relevant atomtypes
	monomer_atomtype = atomtypedefs(0)
	water_atomtype = atomtypedefs(1)
	hydrogen_atomtype = atomtypedefs(2)

	
	DO i = 0, NAtom-2			
		! precompute central atom position for speed
		Posi = Pos(i,:)
		
		Do j = i+1, NAtom-1
			! ignore hydrogens and do typechecking
			if (AtomTypes(i) == hydrogen_atomtype .OR. AtomTypes(j) == hydrogen_atomtype) THEN
				CYCLE
			END IF
			
			IF (AtomTypes(i) == monomer_atomtype .AND. AtomTypes(j) == monomer_atomtype) THEN
				rdf_type = 0 ! monomer-monomer
			ELSE IF (AtomTypes(i) == water_atomtype .AND. AtomTypes(j) == water_atomtype) THEN
				rdf_type = 1 ! water-water
			ELSE IF (AtomTypes(i) == water_atomtype .AND. AtomTypes(j) == monomer_atomtype) THEN
				rdf_type = 2 ! monomer-water
			ELSE IF (AtomTypes(i) == monomer_atomtype .AND. AtomTypes(j) == water_atomtype) THEN
				rdf_type = 2 ! monomer-water
			END IF
			
			
			! Compute pairwise distance
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum-imaging
			rsq = SUM(rij*rij)
			
			! Bin the data
			IF (rsq < Bin_minsq .OR. rsq > Bin_maxsq) THEN
				CYCLE
			END IF
			
			r = sqrt(rsq)
			
			Bin_assignment = AINT((r-Bin_min) * invBin_delta)
			g(Bin_assignment, rdf_type) = g(Bin_assignment, rdf_type) + 1.0	
			
		ENDDO
	ENDDO
	
END SUBROUTINE		



!!! This routine calculates the local density distribution of a central
!!! atomtype with a neighbor atomtype. For atomtypes that are neither, 
!!! the routine returns 0.

SUBROUTINE find_LOCALDENSITY(LD, Pos, NAtom, AtomTypes, Central_atomtype, &
Neigh_atomtype, coeff, LDUpperCuts, NCuts, BoxL)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NCuts, NAtom, Central_atomtype, Neigh_atomtype
	
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1, 0:3) :: coeff
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1) :: LDUpperCuts
	REAL(8), INTENT(IN)  :: BoxL
	
	REAL(8), INTENT(IN), DIMENSION(0:NAtom-1,0:2) :: Pos
	REAL(8), INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
	
	REAL(8), INTENT(INOUT), DIMENSION(0:NAtom-1,0:NCuts-1) :: LD
	!f2py intent(in, out) :: LD
		
	INTEGER :: i,j,k
	REAL(8) :: rsq, phi, c0,c2,c4,c6, R1sq, R2sq, invBoxL
	REAL(8), DIMENSION(0:2) :: rij, Posi
	REAL(8), DIMENSION(0:NCuts-1) :: LDUpperCutsq, LDLowerCutsq
	
	! precompute for speed
	invBoxL = 1.0/BoxL
	DO i = 0, NCuts - 1
		LDUpperCutsq(i) = LDUpperCuts(i) * LDUpperCuts(i)
		LDLowerCutsq(i) = 0.8 * 0.8 * LDUpperCutsq(i)
	ENDDO
	
	DO i = 0, NAtom - 2					
		! precompute central atom position for speed
		Posi = Pos(i,:)
			
		DO j = i+1, NAtom - 1
			! compute pairwise distance
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum imaging
			rsq = SUM(rij*rij)
				
			! calculate local density for all the cutoffs
			DO k = 0, NCuts - 1
				! precompute coefficients and cutoffs for speed
				c0 = Coeff(k,0)
				c2 = Coeff(k,1)
				c4 = Coeff(k,2)
				c6 = Coeff(k,3) 
				R1sq = LDLowerCutsq(k)
				R2sq = LDUpperCutsq(k)
				IF (rsq > R2sq)  THEN
					phi = 0.0
				ELSE IF (rsq <= R1sq) THEN
					phi = 1.0
				ELSE
					phi = c0 + rsq*(c2 + rsq*(c4 + rsq*c6))
				ENDIF
			
				! local density typechecking
				IF (AtomTypes(i) == Central_atomtype .AND. AtomTypes(j) == Neigh_atomtype) THEN
					LD(i,k) = LD(i,k) + phi
				END IF
				
				IF (AtomTypes(i) == Neigh_atomtype .AND. AtomTypes(j) == Central_atomtype) THEN
					LD(j,k) = LD(j,k) + phi 
				END IF
			ENDDO
						
		ENDDO
	ENDDO

END SUBROUTINE	



!!! This routine calculates the # of first shell waters and local density of each 
!!! monomer. It assumes the following atomtypes explicitly: hydrogen - 1, oxygen - 2
!!! monomer - 3. Also it assumes that the position array passed to it has the last N_mon
!!! entries as monomer positions, while the N_water * 3 entries (2 hydrogens and 1 oxygen
!!! per water) come before it. So, the index of the monomer is calculated as :
!!! i - 3*number of waters, where i is the instantaneous loop counter over atoms

SUBROUTINE find_FIRSTSHELLWATERS(ld, fsw, Pos, NAtom, N_water, N_mon, &
AtomTypes, NCuts, LDUpperCuts, FirstShellCut, coeff, BoxL, atomtypedefs)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NAtom, NCuts, N_mon, N_water
	REAL(8), INTENT(IN) :: FirstShellCut
	REAL(8), INTENT(IN), DIMENSION(0:NAtom-1,0:2) :: Pos
	INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
	
	REAL(8), INTENT(IN) :: BoxL
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1,0:3) :: coeff
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1) :: LDUpperCuts
	INTEGER, INTENT(IN), DIMENSION(0:2) :: atomtypedefs
	
	REAL(8), INTENT(INOUT), DIMENSION(0:N_mon - 1, 0:NCuts - 1) :: ld
	REAL(8), INTENT(INOUT), DIMENSION(0:N_mon - 1) :: fsw
	!f2py intent(in,out) :: ld, fsw

	INTEGER :: i,j,k, N_water_total
	REAL(8) :: rsq, phi, FirstShellCutsq, invBoxL
	REAL(8), DIMENSION(0:2) :: Posi, rij 
	REAL(8) :: c0, c2, c4, c6, R1sq, R2sq
	REAL(8), DIMENSION(0:NCuts-1) :: LDUpperCutsq, LDLowerCutsq
	INTEGER :: monomer_atomtype, water_atomtype, hydrogen_atomtype
	
	! Precompute for speed
	N_water_total = 3*N_water
	FirstShellCutsq = FirstShellCut * FirstShellCut
	invBoxL = 1.0/BoxL
	DO i = 0, NCuts-1
		LDUpperCutsq(i) = LDUpperCuts(i) * LDUpperCuts(i)
		LDLowerCutsq(i) = 0.8 * 0.8 * LDUpperCutsq(i) 
	ENDDO
	
	! define relevant atomtypes
	monomer_atomtype = atomtypedefs(0)
	water_atomtype = atomtypedefs(1)
	hydrogen_atomtype = atomtypedefs(2)

		
	DO i = 0, NAtom-1
		! Precompute monomer position for speed
		Posi = Pos(i,:)	
		
		DO j = i+1, NAtom -1
			! ignore hydrogens and oxygen-oxygen cases
			IF (AtomTypes(i) == hydrogen_atomtype .OR. AtomTypes(j) == hydrogen_atomtype) THEN
				CYCLE
			ELSE IF (AtomTypes(i) == water_atomtype .AND. AtomTypes(j) == water_atomtype) THEN
				CYCLE
			END IF
				
			! compute pair-distances
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum-imaging
			rsq = SUM(rij*rij)				
								
			! compute first shell waters
			!! only (i,j) typechecking done since monomer indices are from 5100 to 5125
			IF (AtomTypes(i) == water_atomtype .AND. AtomTypes(j) == monomer_atomtype) THEN
				IF (rsq <= FirstShellCutsq) THEN
					fsw(j-N_water_total) = fsw(j-N_water_total) + 1.00
				END IF			
			END IF
			
			! count monomer local density for each cutoff
			IF (AtomTypes(i) == monomer_atomtype .AND. AtomTypes(j) == monomer_atomtype) THEN
				DO k = 0, NCuts - 1
					! precompute for speed
					R1sq = LDLowerCutsq(k)
					R2sq = LDUpperCutsq(k)
					c0 = coeff(k,0)
					c2 = coeff(k,1)
					c4 = coeff(k,2)
					c6 = coeff(k,3)
				
					IF (rsq >= R2sq) THEN
						phi = 0
					ELSE IF (rsq <= R1sq) THEN
						phi = 1
					ELSE
						phi = c0 + rsq*(c2 + rsq*(c4 + rsq*c6))
					ENDIF
				
					ld(i-N_water_total,k) = ld(i-N_water_total,k) + phi
					ld(j-N_water_total,k) = ld(j-N_water_total,k) + phi
					
				ENDDO
			ENDIF
		
		ENDDO
	ENDDO
	
	
END SUBROUTINE	



!! LD sensitivity tests for pure SPC/E water

SUBROUTINE find_SPCE_FSW(ld, fsw, Pos, NAtom, &
AtomTypes, NCuts, LDUpperCuts, FirstShellCut, coeff, BoxL)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NAtom, NCuts
	REAL(8), INTENT(IN) :: FirstShellCut
	REAL(8), INTENT(IN), DIMENSION(0:NAtom-1,0:2) :: Pos
	INTEGER, INTENT(IN), DIMENSION(0:NAtom-1) :: AtomTypes
	
	REAL(8), INTENT(IN) :: BoxL
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1,0:3) :: coeff
	REAL(8), INTENT(IN), DIMENSION(0:NCuts-1) :: LDUpperCuts
	
	REAL(8), INTENT(INOUT), DIMENSION(0:NAtom - 1, 0:NCuts - 1) :: ld
	REAL(8), INTENT(INOUT), DIMENSION(0:NAtom - 1) :: fsw
	!f2py intent(in,out) :: ld, fsw

	INTEGER :: i,j,k
	REAL(8) :: rsq, phi, FirstShellCutsq, invBoxL
	REAL(8), DIMENSION(0:2) :: Posi, rij 
	REAL(8) :: c0, c2, c4, c6, R1sq, R2sq
	REAL(8), DIMENSION(0:NCuts-1) :: LDUpperCutsq, LDLowerCutsq
	
	FirstShellCutsq = FirstShellCut * FirstShellCut
	invBoxL = 1.0/BoxL
	DO i = 0, NCuts-1
		LDUpperCutsq(i) = LDUpperCuts(i) * LDUpperCuts(i)
		LDLowerCutsq(i) = 0.8 * 0.8 * LDUpperCutsq(i) 
	ENDDO
	
		
	DO i = 0, NAtom-1
		Posi = Pos(i,:)	
		
		DO j = i+1, NAtom -1
			! ignore hydrogens
			IF (AtomTypes(i) == 2 .OR. AtomTypes(j) == 2) THEN
				CYCLE
			END IF
				
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum-imaging
			rsq = SUM(rij*rij)				
								
			IF (AtomTypes(i) == 1 .AND. AtomTypes(j) == 1) THEN
				
				! compute first shell waters
				IF (rsq <= FirstShellCutsq) THEN
					fsw(i) = fsw(i) + 1.00
					fsw(j) = fsw(j) + 1.00
				END IF			
			
				! compute local density for each of the cutoffs
				DO k = 0, NCuts - 1
					R1sq = LDLowerCutsq(k)
					R2sq = LDUpperCutsq(k)
					c0 = coeff(k,0)
					c2 = coeff(k,1)
					c4 = coeff(k,2)
					c6 = coeff(k,3)
				
					IF (rsq >= R2sq) THEN
						phi = 0
					ELSE IF (rsq <= R1sq) THEN
						phi = 1
					ELSE
						phi = c0 + rsq*(c2 + rsq*(c4 + rsq*c6))
					ENDIF
				
					ld(i,k) = ld(i,k) + phi
					ld(j,k) = ld(j,k) + phi
					
				ENDDO
			
			ENDIF
		ENDDO
	ENDDO
		
END SUBROUTINE	









