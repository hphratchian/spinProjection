INCLUDE 'mqc_binary.F03'
INCLUDE 'spinprojection_mod.f03'
      program spinProjection
!
!     This program carries out spin projection on an unrestricted Hartree-Fock
!     (or KS DFT) determinant. The program reads a matrix element file generated
!     from a Gaussian unrestricted HF or DFT calculation.
!
!     This code is based on a previous program written by Hratchian and relies
!     on S2 matrix element code originally written by Lee Thompson.
!
!
!     -H. P. Hratchian, 2023.
!
!
!     USE Connections
!
      use MQC_Binary
      use spinprojection_mod
      use OMP_LIB
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::nCommands,nSwitches,iPrint,nOMP,i,j,k,k1,k2,  &
        nAtoms,nAt3,nBasis,nDetAlpha,nDetBeta,nDetTotal,NDetTT,  &
        nBit_Ints,iAlpha,iBeta,jAlpha,jBeta,nFrozenCore,nFrozenVirtual
      integer(kind=int64),dimension(:),allocatable::atomicNumbers,  &
        tmpVectorInt1,tmpVectorInt2,tmpVectorInt3,tmpVectorInt4
      integer(kind=int64),dimension(:,:),allocatable::permuteAlpha,  &
        permuteBeta
      integer(kind=int64),dimension(:,:,:),allocatable::stringLeftAlpha,  &
        stringLeftBeta,stringRightAlpha,stringRightBeta
      real(kind=real64)::t1,t2,t1A,t2A,Vnn,Escf
      real(kind=real64),dimension(3)::tmp3Vec
      real(kind=real64),dimension(:),allocatable::cartesians,  &
        S2_MatSymm,tmpEVals,tmpVector
      real(kind=real64),dimension(:,:),allocatable::distanceMatrix,  &
        S2_Mat,tmp2NBasisSq,tmpEVecs,tmpMatrix
      character(len=512)::matrixFilename,tmpString
      logical::DEBUG=.false.,gaussianCall,correspondingOrbitals,  &
        ABOverlapTest
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      type(MQC_Variable)::tmpMQCvar
      type(MQC_Scalar)::tmpMQCvar1,tmpMQCvar2,tmpMQCvar3,tmpMQCvar4,  &
        tmpMQCvar5
      type(MQC_Variable)::nEalpha,nEbeta,nEtot,KEnergy,VEnergy,OneElEnergy,  &
        TwoElEnergy,scfEnergy,SSqTmp,SzTotal,MultiplicityTemp
      type(MQC_Variable)::SMatrixAO,SMatrixMOAB,TMatrixAO,VMatrixAO,  &
        HCoreMatrixAO,FMatrixAlpha,FMatrixBeta,PMatrixAlpha,  &
        PMatrixBeta,PMatrixTotal,ERIs,JMatrix,KMatrixAlpha,  &
        CAlpha,CBeta
      type(MQC_R4Tensor)::tmpR4
      Type(MQC_Determinant)::Determinants
!
!     Format Statements
!
 1000 Format(1x,'Enter Program spinProject.')
 1010 Format(3x,'Matrix File: ',A)
 1020 Format(3x,'Print Level: ',I2)
 1030 Format(3x,'Number OMPs: ',I3,/)
 1100 Format(1x,'nAtoms=',I4)
 1200 Format(1x,'Atomic Coordinates (Angstrom)')
 1210 Format(3x,I3,2x,A2,5x,F7.4,3x,F7.4,3x,F7.4)
 1300 Format(1x,'Nuclear Repulsion Energy = ',F20.6)
 2000 Format(1x,60('-'),/,1x,'left  strings: (',I7,' | ',I7,')',  &
        2x,'==> (',B64,' | ',B64,' )',/,  &
        1x,'right strings: (',I7,' | ',I7,')',2x,  &
        '==> (',B64,' | ',B64,' )')
 2010 Format(1x,'S^2 = ',F15.6)
 2100 Format(1x,'Non-Zero Projection (',F10.4,') onto Reference by EVec ',  &
        I10,'. EVal = ',F10.4)
 5000 Format(1x,A,' CPU Time: ',F12.2,' s.')
 8999 Format(/,1x,'END OF TEST PROGRAM scfEnergyTerms.')
!
!
!     Start the job timer and let the user know we are starting the program.
!
      call cpu_time(t1)
      write(IOut,1000)
      call flush(iOut)
!
!     Read the command line arguments.
!
      nCommands = command_argument_count()
      if(nCommands.lt.1)  &
        call mqc_error('No command line arguments provided. At least one command line argument is required giving the input Gaussian matrix element file name.')

!hph+
      write(*,*)
      write(*,*)
      call commandLineArgs_switches(nSwitches,gaussianCall,  &
        correspondingOrbitals,ABOverlapTest,nFrozenCore,nFrozenVirtual,  &
        permuteAlpha,permuteBeta)
      write(*,*)' Hrant - nSwitches = ',nSwitches
      write(*,*)'         gaussianCall          = ',gaussianCall
      write(*,*)'         correspondingOrbitals = ',correspondingOrbitals
      write(*,*)'         nFrozenCore           = ',nFrozenCore
      write(*,*)'         nFrozenVirtual        = ',nFrozenVirtual
      call mqc_print(iOut,permuteAlpha,header='Alpha Permutations')
      call mqc_print(iOut,permuteBeta,header='Beta Permutations')
!hph-

      if(gaussianCall) then
        call commandLineArgs_gaussian(iOut,matrixFilename,iPrint,nOMP)
      else
        call commandLineArgs_direct(nSwitches,matrixFilename,iPrint,nOMP)
      endIf
      write(iOut,*)' matrixFilename: ',TRIM(matrixFilename)
      write(iOut,1020) iPrint
      if(iPrint.eq.-1) then
        iPrint = 10
        DEBUG = .true.
      endIf
      write(iOut,1030) nOMP

!hph+
!      goto 999
!hph-

!
!     Initiate OpenMP.
!
      call omp_set_num_threads(nOMP)
!
!     Open the Gaussian matrix file and load the number of atomic centers.
!
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAtoms = GMatrixFile%getVal('nAtoms')
      write(IOut,1100) nAtoms
!
!     Figure out nAt3, then allocate memory for key arrays.
!
      nAt3 = 3*nAtoms
      Allocate(cartesians(NAt3),atomicNumbers(NAtoms))
!
!     Set the intrinsic nBasis and allocate space for the temp intrinsic nBasis
!     x nBasis matrix.
!
      nBasis = Int(GMatrixFile%getVal('nbasis'))
      Allocate(tmp2NBasisSq(2*nBasis,2*nBasis))
!
!     Load up a few matrices from the matrix file.
!
      call GMatrixFile%getArray('OVERLAP',mqcVarOut=SMatrixAO)
      call GMatrixFile%getArray('KINETIC ENERGY',mqcVarOut=TMatrixAO)
      call GMatrixFile%getArray('CORE HAMILTONIAN ALPHA',mqcVarOut=HCoreMatrixAO)
      call GMatrixFile%getArray('ALPHA FOCK MATRIX',mqcVarOut=FMatrixAlpha)
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA FOCK MATRIX',mqcVarOut=FMatrixBeta)
      else
        FMatrixBeta  = FMatrixAlpha
      endIf
      call GMatrixFile%getArray('ALPHA SCF DENSITY MATRIX',mqcVarOut=PMatrixAlpha)
      call GMatrixFile%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha)
      if(GMatrixFile%isUnrestricted()) then
        call GMatrixFile%getArray('BETA SCF DENSITY MATRIX',mqcVarOut=PMatrixBeta)
        call GMatrixFile%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta)
      else
        PMatrixBeta = PMatrixAlpha
        CBeta = CAlpha
      endIf
      PMatrixTotal = PMatrixAlpha+PMatrixBeta
      VMatrixAO = HCoreMatrixAO-TMatrixAO
!
!     Make permutations for the alpha and beta MO coefficients.
!

      write(*,*)
      write(*,*)
      call CAlpha%print(iOut,header='CAlpha BEFORE Permutations')
      call CBeta%print(iOut,header='CBeta  BEFORE Permutations')
      if(Allocated(permuteAlpha)) then
        do i = 1,SIZE(permuteAlpha,2)
          call MQC_Variable_MatrixPermuteColumns(CAlpha,  &
            permuteAlpha(1,i),permuteAlpha(2,i))
        endDo
      endIf
      if(Allocated(permuteBeta)) then
        do i = 1,SIZE(permuteBeta,2)
          call MQC_Variable_MatrixPermuteColumns(CBeta,  &
            permuteBeta(1,i),permuteBeta(2,i))
        endDo
      endIf
      write(*,*)
      write(*,*)
      call CAlpha%print(iOut,header='CAlpha AFTER  Permutations')
      call CBeta%print(iOut,header='CBeta  AFTER  Permutations')
      write(*,*)
      write(*,*)


!hph      goto 999


!
!     Calculate the number of electrons using <PS>.
!
      nEalpha = Contraction(PMatrixAlpha,SMatrixAO)
      nEbeta  = Contraction(PMatrixBeta,SMatrixAO)
      nEtot   = Contraction(PMatrixTotal,SMatrixAO)
      nEalpha = GMatrixFile%getVal('nalpha')
      nEbeta = GMatrixFile%getVal('nbeta')
      call nEalpha%print(IOut,' <P(Alpha)S>=')
      call nEbeta%print(IOut,' <P(Beta )S>=')
      call nEtot%print(IOut,' <P(Total)S>=')
!
!     Build the total energy.
!
      scfEnergy = MQC(0.5)*  &
        (Contraction(PMatrixAlpha,(HCoreMatrixAO+FMatrixAlpha)) +  &
        Contraction(PMatrixBeta,(HCoreMatrixAO+FMatrixBeta)))
      call scfEnergy%print(iOut,'electronic energy (no Vnn) = ')
!
!     Calculate the 1-electron energy and component pieces of the 1-electron
!     energy. Also, calculate the 2-electron energy.
!
      KEnergy     = Contraction(PMatrixTotal,TMatrixAO)
      VEnergy     = Contraction(PMatrixTotal,VMatrixAO)
      OneElEnergy = Contraction(PMatrixTotal,HCoreMatrixAO)
      TwoElEnergy = scfEnergy-OneElEnergy
      call KEnergy%print(IOut,' <P.K> = ')
      call VEnergy%print(IOut,' <P.V> =')
      call OneElEnergy%print(IOut,' <P.H> = ')
      call TwoElEnergy%print(IOut,' EE, <P.F>-<P.H> = ')
!
!     Load the atommic numbers and Cartesian coordinates into our intrinsic
!     arrays.
!
      atomicNumbers = GMatrixFile%getAtomicNumbers()
      cartesians = GMatrixFile%getAtomCarts()
      cartesians = cartesians*angPBohr
!
!     Print out the atomic numbers and Cartesian coordiantes for each atomic
!     center.
!
      write(IOut,1200)
      do i = 1,NAtoms
        j = 3*(i-1)
        write(IOut,1210) i,mqc_element_symbol(atomicNumbers(i)),  &
          cartesians(j+1),cartesians(j+2),cartesians(j+3)
      endDo
!
!     Form the distance matrix between atomic centers.
!
      Allocate(distanceMatrix(nAtoms,nAtoms))
      do i = 1,nAtoms-1
        distanceMatrix(i,i) = float(0)
        k1 = 3*(i-1)+1
        do j = i+1,NAtoms
          k2 = 3*(j-1)+1
          tmp3Vec = cartesians(k1:k1+2)-cartesians(k2:k2+2)
          distanceMatrix(i,j) = sqrt(dot_product(tmp3Vec,tmp3Vec))
          distanceMatrix(j,i) = distanceMatrix(i,j)
        endDo
      endDo
!
!     Calculate the nuclear-nuclear repulsion energy.
!
      distanceMatrix = distanceMatrix/angPBohr
      Vnn = float(0)
      do i = 1,NAtoms-1
        do j = i+1,NAtoms
          Vnn = Vnn + float(atomicNumbers(i)*atomicNumbers(j))/distanceMatrix(i,j)
        endDo
      endDo
      write(iOut,1300) Vnn
!
!     Put things together and report the SCF energy.
!
      scfEnergy = scfEnergy + MQC(Vnn)
      call scfEnergy%print(IOut,' Total SCF Energy = ')
      call flush(iOut)
!
!     Build the alpha/beta MO overlap matrix and print it out.
!
      if(DEBUG) call SMatrixAO%print(header='AO Overlap Matrix')
      SMatrixMOAB = MatMul(TRANSPOSE(CAlpha),MatMul(SMatrixAO,CAlpha))
      if(DEBUG) call SMatrixMOAB%print(header='CAlpha.S.CAlpha')
      SMatrixMOAB = MatMul(TRANSPOSE(CBeta),MatMul(SMatrixAO,CBeta))
      if(DEBUG) call SMatrixMOAB%print(header='CBeta.S.CBeta')
      SMatrixMOAB = MatMul(TRANSPOSE(CAlpha),MatMul(SMatrixAO,CBeta))
      if(iPrint.ge.0.or.DEBUG)  &
        call SMatrixMOAB%print(header='Alpha-Beta MO Overlap Matrix')
      if(ABOverlapTest) goto 999
!
!     Try building things up for post-SCF like jobs.
!
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - About to call Gen_Det_Str.'
      call flush(iOut)
      tmpMQCvar1 = GMatrixFile%getVal('nbasis')
      tmpMQCvar2 = INT(nEalpha)
      tmpMQCvar3 = INT(nEbeta)
      tmpMQCvar4 = nFrozenCore
      tmpMQCvar5 = nFrozenVirtual
      tmpMQCvar1 = tmpMQCvar1 - tmpMQCvar4 - tmpMQCvar5
      tmpMQCvar2 = tmpMQCvar2 - tmpMQCvar4
      tmpMQCvar3 = tmpMQCvar3 - tmpMQCvar4
      call cpu_time(t1A)
      Call Gen_Det_Str(IOut,2,tmpMQCvar1,tmpMQCvar2,  &
        tmpMQCvar3,Determinants,tmpMQCvar4)
      call cpu_time(t2A)
      write(iOut,5000) 'Gen_Det_Str',t2A-t1A
      nDetAlpha = Bin_Coeff(GMatrixFile%getVal('nbasis')-nFrozenCore-nFrozenVirtual,  &
        INT(nEalpha)-nFrozenCore)
      nDetBeta  = Bin_Coeff(GMatrixFile%getVal('nbasis')-nFrozenCore-nFrozenVirtual,  &
        INT(nEbeta)-nFrozenCore)
      nDetTotal = nDetAlpha*nDetBeta
      write(iOut,*)
      write(iOut,*)' nBasis    = ',INT(GMatrixFile%getVal('nbasis'))
      write(iOut,*)' nEalpha   = ',INT(nEalpha)
      write(iOut,*)' nEbeta    = ',INT(nEbeta)
      write(iOut,*)' nDetAlpha = ',nDetAlpha
      write(iOut,*)' nDetBeta  = ',nDetBeta
      write(iOut,*)' nDetTotal = ',nDetTotal
      write(iOut,*)
      nBit_Ints = ((INT(GMatrixFile%getVal('nbasis'))-nFrozenVirtual)/(Bit_Size(0)-1))+1 
      write(iOut,*)' nBit_Ints = ',nBit_Ints
      write(iOut,*)
      call flush(iOut)

!hph+
!      goto 999
!hph-

      Allocate(S2_Mat(nDetTotal,nDetTotal))
      tmp2NBasisSq(1:NBasis,1:NBasis) = float(0)
      tmp2NBasisSq(NBasis+1:2*NBasis,NBasis+1:2*NBasis) = float(0)
      call MQC_Variable_mqc2intrinsicReal2Array(tmpMatrix,SMatrixMOAB)
      tmp2NBasisSq(1:NBasis,NBasis+1:2*NBasis) = tmpMatrix
      tmp2NBasisSq(NBasis+1:2*NBasis,1:NBasis) = TRANSPOSE(tmpMatrix)
      do i = 1,2*nBasis
        tmp2NBasisSq(i,i) = float(1)
      endDo
      if(iPrint.ge.0.or.DEBUG)  &
        call mqc_print(iOut,tmp2NBasisSq,header='Full MO-MO Overlap')
      call flush(iOut)
!
!     Pre-process the string list combinations for the S2 matrix element
!     formation loops below.
!
      call cpu_time(t1A)
      Allocate(stringLeftAlpha(nDetTotal,nDetTotal,nBit_Ints),  &
        stringLeftBeta(nDetTotal,nDetTotal,nBit_Ints),  &
        stringRightAlpha(nDetTotal,nDetTotal,nBit_Ints),  &
        stringRightBeta(nDetTotal,nDetTotal,nBit_Ints),  &
        tmpVectorInt1(nBit_Ints),tmpVectorInt2(nBit_Ints),  &
        tmpVectorInt3(nBit_Ints),tmpVectorInt4(nBit_Ints))
      call flush(iOut)
      do iBeta = 1,NDetBeta
        tmpVectorInt1 = Determinants%Strings%Beta%vat([iBeta],  &
          [1,nBit_Ints])
        do iAlpha = 1,nDetAlpha
          tmpVectorInt2 = Determinants%Strings%Alpha%vat([iAlpha],  &
            [1,nBit_Ints])
          i = (iBeta-1)*nDetAlpha+iAlpha
          stringLeftBeta(i,1,:) = tmpVectorInt1
          stringLeftAlpha(i,1,:) = tmpVectorInt2
        endDo
      endDo
      call flush(iOut)
      do jBeta = 1,NDetBeta
        tmpVectorInt3 = Determinants%Strings%Beta%vat([jBeta],  &
          [1,nBit_Ints])
        do jAlpha = 1,nDetAlpha
          j = (jBeta-1)*nDetAlpha+jAlpha
          stringRightBeta(1,j,:) = tmpVectorInt3
          tmpVectorInt4 = Determinants%Strings%Alpha%vat([jAlpha],  &
            [1,nBit_Ints])
          stringRightAlpha(1,j,:) = tmpVectorInt4
        endDo
      endDo
      call cpu_time(t2A)
      write(iOut,5000) 'S2 Matrix String Pre-Processing',t2A-t1A
      call flush(iOut)
!
!     Evaluate the S2 matrix elements.
!
      call cpu_time(t1A)
!$OMP PARALLEL DO DEFAULT(NONE),  &
!$OMP SHARED(S2_Mat,stringLeftAlpha,stringLeftBeta,stringRightAlpha,stringRightBeta,  &
!$OMP   tmp2NBasisSq,NDetAlpha,NDetBeta,nBasis,iOut,DEBUG),  &
!$OMP PRIVATE(i,j)
      do i = 1,NDetAlpha*NDetBeta
!        do j = 1,NDetAlpha*NDetBeta
        do j = 1,i
!          S2_Mat(i,j) = S2_Mat_Element_NEW(IOut,2,nBasis,  &
!            stringLeftAlpha(i,1,:),stringLeftBeta(i,1,:),  &
!            stringRightAlpha(1,j,:),stringRightBeta(1,j,:),tmp2NBasisSq,1_int64)
          S2_Mat(i,j) = S2_Mat_Element_NEW(IOut,2,nBasis,  &
            stringLeftAlpha(i,1,:),stringLeftBeta(i,1,:),  &
            stringRightAlpha(1,j,:),stringRightBeta(1,j,:),tmp2NBasisSq)
          if(ABS(S2_Mat(i,j)).lt.(float(1)/float(10000))) S2_Mat(i,j) = float(0)
          S2_Mat(j,i) = S2_Mat(i,j)
!          if(DEBUG) write(iOut,2010) S2_Mat(i,j)
        endDo
      endDo
!$OMP END PARALLEL DO
      call cpu_time(t2A)
      write(iOut,5000) 'S2 Matrix Formation',t2A-t1A
      DeAllocate(stringLeftAlpha,stringLeftBeta,stringRightAlpha,  &
        stringRightBeta)
      if(iPrint.ge.2.or.DEBUG)  &
        call mqc_print(iOut,S2_Mat,header='S2 Matrix (FULL)')
      call flush(iOut)
!
!     Diagonalize S2_Mat and report the eigenvalues.
!
      nDetTT = (nDetTotal*(nDetTotal+1))/2
      Allocate(S2_MatSymm(nDetTT),tmpEVals(nDetTotal),  &
        tmpEVecs(nDetTotal,nDetTotal),tmpVector(3*nDetTotal))
      k = 0
      do i = 1,nDetTotal
        do j = i,nDetTotal
          k = k+1
          S2_MatSymm(k) = S2_Mat(i,j)
        endDo
      endDo
      if(DEBUG) call mqc_print(iOut,mqc_matrixSymm2Full(S2_MatSymm,'L'),header='S2 Matrix from Symm')
      write(iOut,*)' k = ',k
      write(iOut,*)' nDetTT = ',nDetTT
      write(iOut,*)
      write(iOut,*)
!hph      goto 999
      write(iOut,*)' Calling DSPEV.'
      call cpu_time(t1A)
      Call DSPEV('V','L',nDetTotal,S2_MatSymm,tmpEVals,tmpEVecs,  &
        NDetTotal,tmpVector,i)
      call cpu_time(t2A)
      write(iOut,5000) 'S2 Diagonalization',t2A-t1A
      write(iOut,*)' After DSPEV, i = ',i
      call flush(iOut)
      do i = 1,nDetTotal
        if(ABS(tmpEVals(i)).lt.(float(1)/float(10000))) tmpEVals(i) = float(0)
      endDo
      call mqc_print(iOut,tmpEVals,header='S2 EVals')
      if(iPrint.ge.2.or.DEBUG)  &
        call mqc_print(iOut,tmpEVecs,header='S2 EVecs')
      if(DEBUG) call mqc_print(iOut,MatMul(Transpose(tmpEVecs),MatMul(S2_Mat,tmpEVecs)),header='EVecs.S2.EVecs')
      if(DEBUG) call mqc_print(iOut,MatMul(Transpose(tmpEVecs),tmpEVecs),header='EVecs.EVecs')
      write(iOut,*)
      write(iOut,*)' nDetTotal = ',nDetTotal
      call flush(iOut)
      do i = 1,NDetTotal
        if(ABS(tmpEVecs(nDetTotal,i)).gt.(float(5)/float(100))) then
          write(iOut,2100) tmpEVecs(nDetTotal,i),i,tmpEvals(i)
        endIf
      endDo
!
!     Try to figure out the range of Sz values found in the S^2 eigenvalue
!     range.
!
      SzTotal = MQC(0.5)*(nEAlpha-nEBeta)
      MultiplicityTemp = nEAlpha-nEBeta+MQC(1)
      call SzTotal%print(header='Sz=')
      call MultiplicityTemp%print(header='Multiplicity=')
      call sumOverSpinStates(iOut,INT(MultiplicityTemp),tmpEVals,  &
        tmpEVecs(nDetTotal,:))
!
  999 Continue
      write(iOut,8999)
      call cpu_time(t2)
      write(iOut,5000) 'TOTAL JOB',t2-t1
      end program spinProjection
