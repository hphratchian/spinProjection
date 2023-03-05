      module spinprojection_mod
!
!     This module supports the program scfEnergyTerms.
!
!     -H. P. Hratchian, 2020, 2021.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
!
!
!     Module Procedures
!
      CONTAINS
!
!
      subroutine commandLineArgs_switches(nSwitches,gaussianCall,  &
        correspondingOrbitals,ABOverlapTest,nFrozenCore,nFrozenVirtual,  &
        permuteAlpha,permuteBeta)
!
!     This routine reads the command line arguments for switches, which are all
!     required BEFORE the standard command line arguments. This routine returns
!     the number of command line arguments that are switches (or related to
!     switches), which can be used by the calling program to offset the
!     positions of non-switch command line arguments.
!
!
      implicit none
      integer(kind=int64),intent(out)::nSwitches,nFrozenCore,  &
        nFrozenVirtual
      integer(kind=int64),dimension(:,:),allocatable,intent(out)::permuteAlpha,  &
        permuteBeta
      logical,intent(out)::gaussianCall,correspondingOrbitals,ABOverlapTest
      integer(kind=int64)::i,nCommands,iPermute,nPermute
      character(len=512)::tmpString
      logical::foundSwitch
!
!     Default the output arguments.
!
      gaussianCall = .False.
      correspondingOrbitals = .False.
      ABOverlapTest = .False.
      nFrozenCore = 0
      nFrozenVirtual = 0
!
!     Do the work...
!
      nCommands = command_argument_count()
      i = 1
      foundSwitch = .True.
      do while(i.lt.nCommands)
        call get_command_argument(i,tmpString)
        if(tmpString(1:1).ne.'-') foundSwitch = .False.
        if(.not.foundSwitch) exit
        write(*,*)' Hrant - Found Switch: ',TRIM(tmpString)
        select case(TRIM(tmpString))
        case('-g')
          gaussianCall = .true.
          i = i+1
        case('-c')
          correspondingOrbitals = .True.
          i = i+1
        case('-overlap')
          ABOverlapTest = .True.
          i = i+1
        case('-fc')
          call mqc_get_command_argument_integer(i+1,nFrozenCore)
          i = i+2
        case('-fv')
          call mqc_get_command_argument_integer(i+1,nFrozenVirtual)
          i = i+2
        case('-p')
          call mqc_get_command_argument_integer(i+1,nPermute)
          i = i+2
          Allocate(permuteAlpha(2,nPermute),permuteBeta(2,nPermute))
          do iPermute=1,nPermute
            call mqc_get_command_argument_integer(i,permuteAlpha(1,iPermute))
            call mqc_get_command_argument_integer(i+1,permuteAlpha(2,iPermute))
            call mqc_get_command_argument_integer(i,permuteBeta(1,iPermute))
            call mqc_get_command_argument_integer(i+1,permuteBeta(2,iPermute))
            i = i+2
          endDo
        case('-pa')
          call mqc_get_command_argument_integer(i+1,nPermute)
          i = i+2
          Allocate(permuteAlpha(2,nPermute))
          do iPermute=1,nPermute
            call mqc_get_command_argument_integer(i,permuteAlpha(1,iPermute))
            call mqc_get_command_argument_integer(i+1,permuteAlpha(2,iPermute))
            i = i+2
          endDo
        case('-pb')
          call mqc_get_command_argument_integer(i+1,nPermute)
          i = i+2
          Allocate(permuteBeta(2,nPermute))
          do iPermute=1,nPermute
            call mqc_get_command_argument_integer(i,permuteBeta(1,iPermute))
            call mqc_get_command_argument_integer(i+1,permuteBeta(2,iPermute))
            i = i+2
          endDo
        case default
          call mqc_error('Unknown Command Line Switch Found.')
!hph          i = i+1
        end select
      end do
      nSwitches = i-1
!
      return
      end subroutine commandLineArgs_switches
!
!
      subroutine commandLineArgs_direct(nSwitches,matrixFilename,iPrint,  &
        nOMP)
!
!     This routine reads the command line and returns the name of the matrix
!     element file, the print level flag, and the requested number of openMP
!     processes. This routine is meant for use in direct runs by a user.
!
!     The input dummy argument nSwitches is used to offset the command line
!     arguments where the matrix filename, print flag, and number of OpenMP
!     threads are expected.
!
!
      implicit none
      integer(kind=int64),intent(in)::nSwitches
      character(len=512),intent(out)::matrixFilename
      integer(kind=int64),intent(out)::iPrint,nOMP
      integer(kind=int64)::nCommands
      character(len=512)::tmpString
!
!     Do the work...
!
      nCommands = command_argument_count()-nSwitches
      call get_command_argument(1+nSwitches,matrixFilename)
      iPrint = 0
      nOMP = 1
      if(nCommands.ge.2) then
        call get_command_argument(2+nSwitches,tmpString)
        read(tmpString,*) iPrint
      endIf
      if(nCommands.ge.3) then
        call get_command_argument(3+nSwitches,tmpString)
        read(tmpString,*) nOMP
        if(nOMP.le.0) call mqc_error('OMP number must be >= 1.')
      endIf
      if(nCommands.gt.3)  &
        call mqc_error('More than 3 command line arguments provided.')
!
      return
      end subroutine commandLineArgs_direct
!
!
      subroutine commandLineArgs_gaussian(IOut,matrixFilename,iPrint,nOMP)
!
!     This routine reads the command line and returns the name of the matrix
!     element file, the print level flag, and the requested number of openMP
!     processes in cases when the program is called directly from Gaussian.
!     Additionally, this subroutine connects file unit iOut to the Gaussian
!     message file so that the output from this program will be put into the
!     Gaussian log file.
!
!
      implicit none
      integer(kind=int64),intent(in)::iOut
      character(len=512),intent(out)::matrixFilename
      integer(kind=int64),intent(out)::iPrint,nOMP
      integer(kind=int64)::nCommands,iOff,iError
      character(len=512)::tmpString
!
!     Do the work...
!
      iOff = 1
      nCommands = command_argument_count()
      call get_command_argument(4+iOff,tmpString)
      Open(Unit=iOut,File=TRIM(tmpString),Status='unknown',IOStat=IError)
      call get_command_argument(6+iOff,matrixFilename)
      iPrint = 0
      nOMP = 1
!      if(nCommands.ge.2) then
!        call get_command_argument(2,tmpString)
!        read(tmpString,*) iPrint
!      endIf
!      if(nCommands.ge.3) then
!        call get_command_argument(3,tmpString)
!        read(tmpString,*) nOMP
!        if(nOMP.le.0) call mqc_error('OMP number must be >= 1.')
!      endIf
!      if(nCommands.gt.3)  &
!        call mqc_error('More than 3 command line arguments provided.')
!
      return
      end subroutine commandLineArgs_gaussian
!
!PROCEDURE sumOverSpinStates
      subroutine sumOverSpinStates(iOut,lowestMultiplicity,SSqList,  &
        populationOverDets)
!
!     This subroutine takes the lowest expected multiplicity for a system and
!     the list of S^2 eigenvalues, and determines how many states of each spin
!     multiplicity are present in the list.
!
!
      implicit none
      integer(kind=int64),intent(in)::iOut,lowestMultiplicity
      real(kind=real64),dimension(:),intent(in)::SSqList,populationOverDets
      integer(kind=int64)::i,j,nStates,nMultips,multiplicityLow,  &
        multiplicityHigh,multiplicityCurrent
      integer(kind=real64),dimension(:,:),allocatable::multiplicityList
      real(kind=real64)::realTmp
      real(kind=real64),dimension(:),allocatable::populationsOverMultiplicities
!
 1000 format(/,1x,'Number of Multiplicities: ',I5,/,  &
        3x,'Lowest  Multiplicity in S**2: ',i5,/,  &
        3x,'Highest Multiplicity in S**2: ',i5)
 1100 format(1x,i3,3x,f10.6,3x,f10.6,3x,i3)
!
      nStates = SIZE(SSqList)
!
!     Figure out the lowest and highest multiplicities in SSqList. SSqList is
!     assumed to be ordered from lowest to highest value.
!
      multiplicityLow  = INT(SQRT(1.0+4.0*SSqList(1))+0.5)
      multiplicityHigh = INT(SQRT(1.0+4.0*SSqList(nStates))+0.5)
      nMultips = 1 + (multiplicityHigh - multiplicityLow)/2
      write(iOut,1000) nMultips,multiplicityLow,multiplicityHigh
      Allocate(multiplicityList(2,nMultips),  &
        populationsOverMultiplicities(nMultips))
      multiplicityCurrent = multiplicityLow
      do i = 1,nMultips
        multiplicityList(1,i) = multiplicityCurrent
        multiplicityList(2,i) = 0
        multiplicityCurrent = multiplicityCurrent+2
      endDo
      call mqc_print(iOut,multiplicityList,header='Multiplicity List (Initial)')
!
!
!     Loop over SSqList elements and evaluate 1+4*S^2 to see if they're integer
!     values.
      populationsOverMultiplicities = 0.0
      do i = 1,SIZE(SSqList)
        multiplicityCurrent = INT(SQRT(1.0+4.0*SSqList(i))+0.5)
        j = (multiplicityCurrent-multiplicityLow)/2 + 1
        multiplicityList(2,j) = multiplicityList(2,j) + 1
        populationsOverMultiplicities(j) = populationsOverMultiplicities(j)   &
          + populationOverDets(i)**2
!hph        write(iOut,1100) i,SSqList(i),4.0*SSqList(i),multiplicityCurrent
      endDo
      call mqc_print(iOut,multiplicityList,header='Multiplicity List (Final)')
      call mqc_print(iOut,populationsOverMultiplicities,header='Sum of Squares of Amplitudes Over Spin Subspaces')
      write(iOut,*)' Hrant - SUM = ',SUM(populationsOverMultiplicities)
!
      return
      end subroutine sumOverSpinStates
!
!
      subroutine formFock(nBasis,density,ERIs,coulomb)
!
!     This subroutine forms a Coulomb matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::coulomb
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempCoulomb
!
      call ERIs%print(IOut,' In formCoulomb: ERIs=',blankAtTop=.True.)
!
!     Do the work...
!
      tempCoulomb = float(0)
      do iSigma    = 1,nBasis
        do iLambda = 1,nBasis
          do iNu   = 1,nBasis
            do iMu = 1,nBasis
              write(IOut,'(/,1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3,2x,I3,2x,I3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])),  &
                float(density%getVal([iLambda,iSigma])),  &
                MQC_Variable_getArrayPosition(ERIs,[iMu,iNu,iLambda,iSigma])
              tempCoulomb(iMu,iNu) = tempCoulomb(iMu,iNu) +  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      coulomb = tempCoulomb
!
      return
      end subroutine formFock
!
!
      subroutine formCoulomb(nBasis,density,ERIs,coulomb,initialize)
!
!     This subroutine forms a Coulomb matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::coulomb
      logical,optional::initialize
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempCoulomb
      logical::init
!
      call ERIs%print(IOut,' In formCoulomb: ERIs=',blankAtTop=.True.)
!
!     Do initial set-up work.
!
      init = .true.
      if(Present(initialize)) init = initialize
!
!     Work through the integral loops to build Coulomb contributions.
!
      if(init) tempCoulomb = float(0)
      do iSigma    = 1,nBasis
        do iLambda = 1,nBasis
          do iNu   = 1,nBasis
            do iMu = 1,nBasis
              write(IOut,'(/,1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3,2x,I3,2x,I3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])),  &
                float(density%getVal([iLambda,iSigma])),  &
                MQC_Variable_getArrayPosition(ERIs,[iMu,iNu,iLambda,iSigma])
              tempCoulomb(iMu,iNu) = tempCoulomb(iMu,iNu) +  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      coulomb = tempCoulomb
!
      return
      end subroutine formCoulomb
!
!
      subroutine formExchange(nBasis,density,ERIs,exchange,initialize)
!
!     This subroutine forms an Exchange matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::exchange
      logical,optional::initialize
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempExchange
      logical::init
!
!     Do initial set-up work.
!
      init = .true.
      if(Present(initialize)) init = initialize
!
!     Work through the integral loops to build Exchange contributions.
!
      if(init) tempExchange = float(0)
      do iMu = 1,nBasis
        do iNu = 1,nBasis
          do iLambda = 1,nBasis
            do iSigma = 1,nBasis
              write(IOut,'(1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iSigma,iLambda,iNu])),  &
                float(density%getVal([iLambda,iSigma]))
              tempExchange(iMu,iNu) = tempExchange(iMu,iNu) -  &
                float(ERIs%getVal([iMu,iSigma,iLambda,iNu])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      exchange = tempExchange
!
      return
      end subroutine formExchange

!
!=====================================================================
!
!     PROCEDURE S2_Mat_Element_NEW
! 
      Function S2_Mat_Element_NEW(IOut,IPrint,NBasis,Alpha_String_1,Beta_String_1, &
        Alpha_String_2,Beta_String_2,MO_Overlap,nBit_IntsIn)
!
!     This function returns the CI S**2 matrix elements used for computing S**2
!     values of CI vectors for a given alpha and beta string combination.
!     The MO overlap matrix is required
!
!     Variable Declarations...
!
      Implicit None
      Integer(kind=int64),optional,intent(in)::NBit_IntsIn
      Integer(kind=int64)::IOut,IPrint,NBasis,IPos,JPos,IDiff,Det_Diff,NAlpha,NBeta, &
        IOcc,JOcc,KOcc,LOcc,Mat_Sign,Alpha_Diff_Cnt,Beta_Diff_Cnt,NBit_Ints, &
        I,J,II,JJ
      real(kind=real64)::S2_Mat_Element_NEW,Zero=0.0d0,Quarter=0.25d0,ABTerm,One=1.0d0
      Integer(kind=int64),Dimension(4)::Orbs,Spin,Det
      Integer(kind=int64),Dimension(:)::Alpha_String_1,Alpha_String_2,Beta_String_1, &
        Beta_String_2
      Integer(kind=int64),Dimension(:),Allocatable::Alpha_Diff,Beta_Diff
      Real(kind=real64),Dimension(:,:),Allocatable::MO_Overlap
!
!      Write(IOut,*) 'alpha 1:'
!      Write(IOut,'(B64)') Alpha_String_1
!      Write(IOut,*) 'alpha 2:'
!      Write(IOut,'(B64)') Alpha_String_2
!      Write(IOut,*) 'alpha XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Alpha_String_1,Alpha_String_2)   
!      Write(IOut,*) 'NDifa=', PopCnt(IEOR(Alpha_String_1,Alpha_String_2))
!      Write(IOut,*)
!      Write(IOut,*) 'beta 1:'
!      Write(IOut,'(B64)') Beta_String_1
!      Write(IOut,*) 'beta 2:'
!      Write(IOut,'(B64)') Beta_String_2
!      Write(IOut,*) 'beta XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Beta_String_1,Beta_String_2)   
!      Write(IOut,*) 'NDifb=', PopCnt(IEOR(Beta_String_1,Beta_String_2))
!      Write(IOut,*)
!
      if(PRESENT(NBit_IntsIn)) then
        NBit_Ints = NBit_IntsIn
      else
        NBit_Ints = (NBasis/Bit_Size(0))+1 
      endIf
      Allocate(Alpha_Diff(NBit_Ints),Beta_Diff(NBit_Ints))
      Det_Diff = 0
      Alpha_Diff_Cnt = 0
      Beta_Diff_Cnt = 0
      NAlpha = 0
      NBeta = 0
      Do I = 1,NBit_Ints
        Alpha_Diff(I) = IEOR(Alpha_String_1(I),Alpha_String_2(I))
!        Write(IOut,*) 'Alpha Diff',I,':'
!        Write(IOut,'(B64)') Alpha_Diff(I)
!        Write(IOut,*) '-------------'
        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(I)) 
        NAlpha = NAlpha + PopCnt(Alpha_String_1(I))
        Beta_Diff(I) = IEOR(Beta_String_1(I),Beta_String_2(I))
!        Write(IOut,*) 'Beta Diff',I,':'
!        Write(IOut,'(B64)') Beta_Diff(I)
!        Write(IOut,*) '-------------'
        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(I))
        NBeta = NBeta + PopCnt(Beta_String_1(I))
      EndDo
!      Write(IOut,*)'Alpha_Diff_Cnt:',Alpha_Diff_Cnt,'Beta_Diff_Cnt:',Beta_Diff_Cnt
      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2

      If(Mod(Alpha_Diff_Cnt,2).ne.0.or.Mod(Beta_Diff_Cnt,2).ne.0) then
        Write(IOut,*) "ERROR: S2_Mat_Element_NEW has been handed spin non-conserving &
        determinants"
        Call Exit()
      EndIf
      Select Case (Det_Diff)
!
        Case(3:)
          S2_Mat_Element_NEW = Zero 
          Return
!
        Case(2)
          IDiff = 1
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Call Print_Vector(IOut,Orbs,'Orbs')
!          Call Print_Vector(IOut,Spin,'Spin')
!          Call Print_Vector(IOut,Det,'Det')
!
          IOcc = 0
          Do IPos = Orbs(1)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(1).eq.0) then
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            ElseIf(Spin(1).eq.1) then
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'IOcc:',IOcc
          JOcc = 0
          Do IPos = Orbs(2)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(2).eq.0) then
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            ElseIf(Spin(2).eq.1) then
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'JOcc:',JOcc
          KOcc = 0
          Do IPos = Orbs(3)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(3).eq.0) then
              If(Det(3).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) KOcc = KOcc + 1
              EndIf
            ElseIf(Spin(3).eq.1) then
              If(Det(3).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) KOcc = KOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'KOcc:',KOcc
          LOcc = 0
          Do IPos = Orbs(4)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(4).eq.0) then
              If(Det(4).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) LOcc = LOcc + 1
              EndIf
            ElseIf(Spin(4).eq.1) then
              If(Det(4).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) LOcc = LOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'LOcc:',LOcc
!          Mat_Sign = -1
!          Mat_Sign = (-1)**(IOcc+JOcc+KOcc+LOcc-3)
!          Write(IOut,*) 'Permutations:',(IOcc+JOcc+KOcc+LOcc-3)
          Mat_Sign = (-1)**(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Permutations:',(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Mat_Sign:',Mat_Sign
!
          If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
            If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(2),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Element_NEW = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(2),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Element_NEW = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Element_NEW = Zero
            EndIf
          ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(3),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Element_NEW = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(3),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Element_NEW = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Element_NEW = Zero
            EndIf
          ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(4),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Element_NEW = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Element_NEW = Mat_Sign*MO_Overlap(Orbs(4),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Element_NEW = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Element_NEW = Zero
            EndIf
          EndIf
!          Write(IOut,*) 'S2_Mat_Element_NEW:',S2_Mat_Element_NEW

          Return

        Case(1)
          IDiff = 1
!          Allocate(Orbs(2),Spin(2))
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!          Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!          Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!
          S2_Mat_Element_NEW = Zero 
          If(Spin(1).ne.Spin(2)) then
            S2_Mat_Element_NEW = Zero
!
          ElseIf(Spin(1).eq.0) then
!
            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Beta_String_1(I),J).eq..True.) then
                S2_Mat_Element_NEW = S2_Mat_Element_NEW + MO_Overlap(IPos+1+NBasis,Orbs(1)) * MO_Overlap(Orbs(2),IPos+1+NBasis)
              EndIf
            EndDo
!            S2_Mat_Element_NEW = - S2_Mat_Element_NEW
!            Write(IOut,*) 'Permutations:',(2*NAlpha+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NAlpha+1-IOcc-JOcc)
            S2_Mat_Element_NEW = (-1)**(2*NAlpha+1-IOcc-JOcc) * S2_Mat_Element_NEW
!            S2_Mat_Element_NEW = (-1)**(IOcc+JOcc-1) * S2_Mat_Element_NEW
!
          ElseIf(Spin(1).eq.1) then

            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                S2_Mat_Element_NEW = S2_Mat_Element_NEW + MO_Overlap(IPos+1,Orbs(1)+NBasis) * MO_Overlap(Orbs(2)+NBasis,IPos+1)
              EndIf
            EndDo
!            S2_Mat_Element_NEW = - S2_Mat_Element_NEW
!            Write(IOut,*) 'Permutations:',(2*NBeta+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NBeta+1-IOcc-JOcc)
            S2_Mat_Element_NEW = (-1)**(2*NBeta+1-IOcc-JOcc) * S2_Mat_Element_NEW
!            S2_Mat_Element_NEW = (-1)**(IOcc+JOcc-1) * S2_Mat_Element_NEW

          EndIf

!          Write(IOut,*) 'S2_Mat_Element_NEW:',S2_Mat_Element_NEW

          Return
!
        Case(0)
          ABTerm = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = 0, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0)) 
              If((BTest(Alpha_String_2(I),J).eq..True.).and.(BTest(Beta_String_2(II),JJ).eq..True.)) then
                ABTerm = ABTerm + MO_Overlap(IPos+1,JPos+1+NBasis)*MO_Overlap(JPos+1+NBasis,IPos+1) 
              EndIf
            EndDo
          EndDo
          S2_Mat_Element_NEW = Quarter*((NAlpha-NBeta)**2+2*(NAlpha+NBeta)) - ABTerm
!          Write(IOut,*) 'S2_Mat_Element_NEW:',S2_Mat_Element_NEW
!
          Return
!
      End Select
!
      End Function S2_Mat_Element_NEW  

!
!
      end module spinprojection_mod
