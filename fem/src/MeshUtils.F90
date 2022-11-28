!*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh manipulation utilities for *Solver - routines
!------------------------------------------------------------------------------

MODULE MeshUtils

    USE LoadMod
    USE MeshBasics
    USE MortarUtils
    USE ElementUtils
    USE ElementDescription
    USE Interpolation
    USE ParallelUtils
    USE Types
    USE InterpVarToVar
    USE MatrixAssembly, ONLY : mGetElementDofs, mGetBoundaryIndexesFromParent
    IMPLICIT NONE

CONTAINS



!------------------------------------------------------------------------------
! This version of creating def_dofs arrays has limited abilities since it does not
! support element family flags (cf. the subroutine GetDefs in ModelDescription). 
! There is no need for calling this unless the element definition is given in an 
! equation section or a matc function is used to evaluate the order of p-basis,
! since otherwise the subroutine GetDefs has done the necessary work.
! TO DO: Have just one subroutine for writing def_dofs arrays ?
!------------------------------------------------------------------------------
   SUBROUTINE GetMaxDefs(Model, Mesh, Element, ElementDef, SolverId, BodyId, Def_Dofs)
!------------------------------------------------------------------------------
     CHARACTER(*) :: ElementDef
     TYPE(Model_t) :: Model
     TYPE(MEsh_t) :: Mesh
     TYPE(Element_t) :: Element
     INTEGER :: SolverId, BodyId, Def_Dofs(:,:)

     TYPE(ValueList_t), POINTER :: Params
     INTEGER :: i, j,k,l, n, slen, Family
     INTEGER, POINTER :: Body_Dofs(:,:)
     LOGICAL  :: stat, Found
     REAL(KIND=dp) :: x,y,z
     TYPE(Solver_t), POINTER  :: Solver
     CHARACTER(MAX_NAME_LEN) :: str, RESULT

     TYPE(ValueList_t), POINTER :: BodyParams
     CHARACTER(MAX_NAME_LEN) :: ElementDefBody
     
     BodyParams => Model % Bodies(BodyId) % Values

     ElementDefBody=ListGetString(BodyParams,'Solver '//TRIM(i2s(SolverId))//': Element',Found )
     IF (Found) THEN
       CALL Info('GetMaxDefs','Element found for body '//TRIM(i2s(BodyId))//' with solver '//TRIM(i2s(SolverId)), Level=5) 
       CALL Info('GetMaxDefs','Default element type is: '//ElementDef, Level=5)
       CALL Info('GetMaxDefs','New element type for this body is now: '//ElementDefBody, Level=5)
       ElementDef=ElementDefBody
     END IF

     Solver => Model % Solvers(SolverId)
     Params => Solver % Values

     IF ( .NOT. ALLOCATED(Solver % Def_Dofs) ) THEN
       ALLOCATE(Solver % Def_Dofs(10,Model % NumberOfBodies,6))
       Solver % Def_Dofs=-1
       Solver % Def_Dofs(:,:,1)=1
     END IF
     Body_Dofs => Solver % Def_Dofs(1:8,BodyId,:)

     j = INDEX(ElementDef, '-') ! FIX this to include elementtypewise defs...
     IF ( j>0 ) THEN
       CALL Warn('GetMaxDefs', &
           'Element set flags not supported, move element definition to a solver section')
       RETURN
     END IF

     j = INDEX( ElementDef, 'n:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l
       Body_Dofs(:,1) = l
       Def_Dofs(:,1) = MAX(Def_Dofs(:,1), l)
     END IF
          
      j = INDEX( ElementDef, 'e:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,2) = l
        Def_Dofs(1:8,2) = MAX(Def_Dofs(1:8,2), l )
      END IF
          
      j = INDEX( ElementDef, 'f:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,3) = l
        Def_Dofs(1:8,3) = MAX(Def_Dofs(1:8,3), l )
      END IF
          
      j = INDEX( ElementDef, 'd:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l

        ! Zero value triggers discontinuous approximation,
        ! substitute the default negative initialization value to avoid troubles:
        IF (l == 0) l = -1

        Body_Dofs(:,4) = l
        Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4), l )
      ELSE 
        IF ( ListGetLogical( Solver % Values, &
            'Discontinuous Galerkin', stat ) ) THEN
          Body_Dofs(:,4) = 0
          Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4),0 )
        END IF
      END IF
          
      j = INDEX( ElementDef, 'b:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(1:8,5) = l
        Def_Dofs(1:8,5) = MAX(Def_Dofs(1:8,5), l )
      END IF
          
      j = INDEX( ElementDef, 'p:' )
      IF ( j>0 ) THEN
        IF ( ElementDef(j+2:j+2) == '%' ) THEN
          n = Element % TYPE % NumberOfNodes
          x = SUM(Mesh % Nodes % x(Element % NodeIndexes))/n
          y = SUM(Mesh % Nodes % y(Element % NodeIndexes))/n
          z = SUM(Mesh % Nodes % z(Element % NodeIndexes))/n
!          WRITE( str, * ) 'cx= ',TRIM(i2s(Element % ElementIndex)),x,y,z
          WRITE( str, * ) 'cx= ',TRIM(i2s(Element % BodyId)),x,y,z
          str = TRIM(str) // '; ' // TRIM(ElementDef(j+3:))//'(cx)'
          slen = LEN_TRIM(str)
          CALL matc(str,RESULT,slen)
          READ(RESULT(1:slen),*) x

          Def_Dofs(1:8,6)  = MAX(Def_Dofs(1:8,6),NINT(x))
          Family = Element % TYPE % ElementCode / 100
          Body_Dofs(Family, 6) = &
              MAX(Body_Dofs(Family, 6), NINT(x))
        ELSE
          READ( ElementDef(j+2:), * ) l
          Body_Dofs(:,6) = l
          Def_Dofs(1:8,6) = MAX(Def_Dofs(1:8,6), l )
        END IF
      END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetMaxDefs
!------------------------------------------------------------------------------


!> Create a discontinuous mesh over requested boundaries.
!> The nodes are duplicated in order to facilitate the discontinuity.
!> The duplicate nodes are not created by default if the connectivity 
!> of the nodes is needed by other bulk elements than those directly 
!> associated with the discontinuous boundaries. 
!------------------------------------------------------------------------------
 SUBROUTINE CreateDiscontMesh( Model, Mesh, DoAlways )

   TYPE(Model_t) :: Model
   TYPE(Mesh_t), POINTER :: Mesh
   LOGICAL, OPTIONAL :: DoAlways

   INTEGER, POINTER :: DisContPerm(:)
   LOGICAL, ALLOCATABLE :: DisContNode(:), DisContElem(:), ParentUsed(:), &
       MovingNode(:), StayingNode(:)
   LOGICAL :: Found, DisCont, GreedyBulk, GreedyBC, Debug, DoubleBC, UseTargetBodies, &
       UseConsistantBody, LeftHit, RightHit, Moving, Moving2, Set, Parallel
   INTEGER :: i,j,k,l,n,m,t,bc
   INTEGER :: NoNodes, NoDisContElems, NoDisContNodes, &
       NoBulkElems, NoBoundElems, NoParentElems, NoMissingElems, &
       DisContTarget, NoMoving, NoStaying, NoStayingElems, NoMovingElems, &
       NoUndecided, PrevUndecided, NoEdges, Iter, ElemFamily, DecideLimit, &
       ActiveBCs, CandA, CandB, RightBody, LeftBody, ConflictElems
   INTEGER, TARGET :: TargetBody(1)
   INTEGER, POINTER :: Indexes(:),ParentIndexes(:),TargetBodies(:)
   TYPE(Element_t), POINTER :: Element, LeftElem, RightElem, ParentElem, OtherElem
   CHARACTER(MAX_NAME_LEN) :: DiscontFlag
   LOGICAL :: CheckForHalo
   LOGICAL, POINTER :: HaloNode(:)
   TYPE(ValueList_t), POINTER :: BCList
   LOGICAL :: DoneThisAlready = .FALSE.
   CHARACTER(*), PARAMETER :: Caller = 'CreateDiscontMesh'

   IF(.NOT.PRESENT(DoAlways)) THEN
     IF (DoneThisAlready) RETURN
   ELSE 
     IF(.NOT.DoAlways) THEN
       IF (DoneThisAlready) RETURN
     END IF
   END IF
   DoneThisAlready = .TRUE.

   Discont = .FALSE.
   DoubleBC = .FALSE.
   ActiveBCs = 0
   DO bc = 1,Model % NumberOfBCs
     DisCont = ListGetLogical( Model % BCs(bc) % Values,'Discontinuous Boundary',Found )
     ! If the target boundary / periodic bc / mortar bc is zero
     ! it refers to itself. Otherwise the boundary will be doubled.
     IF( DisCont ) THEN
       i = ListGetInteger( Model % BCs(bc) % Values,'Discontinuous BC',Found )
       j = ListGetInteger( Model % BCs(bc) % Values,'Periodic BC',Found )
       k = ListGetInteger( Model % BCs(bc) % Values,'Mortar BC',Found )
       l = ListGetInteger( Model % BCs(bc) % Values,'Contact BC',Found )
       DoubleBC = ( i + j + k + l > 0 )
       ActiveBCs = ActiveBCs + 1
       BCList => Model % BCs(bc) % Values
     END IF
   END DO
   IF(ActiveBCs == 0 ) RETURN
   
   CALL Info(Caller,'Creating discontinuous boundaries')

   IF( ActiveBCs > 1 ) THEN
     CALL Warn(Caller,'Be careful when using more than one > Discontinuous Boundary < !')
   END IF

   Parallel = ( ParEnv % PEs > 1 )

   NoNodes = Mesh % NumberOfNodes
   NoBulkElems = Mesh % NumberOfBulkElements
   NoBoundElems = Mesh % NumberOfBoundaryElements
   
   ALLOCATE( DisContNode(NoNodes))
   ALLOCATE( DisContElem(NoBoundElems))
   ALLOCATE( ParentUsed(NoBulkElems))
   DisContNode = .FALSE.
   DisContElem = .FALSE.
   ParentUsed = .FALSE.
   NoDisContElems = 0
   NoMissingElems = 0


   ! Check whether we need to skip some elements and nodes on the halo boundary 
   ! We might not want to create additional nodes on the nodes that are on the halo only 
   ! since they just would create further need for new halo...
   CheckForHalo = ListGetLogical( Model % Simulation,'No Discontinuous Halo',Found ) 
   IF(.NOT. Found ) CheckForHalo = .TRUE.
   IF( CheckForHalo ) THEN
     HaloNode => NULL()
     CALL MarkHaloNodes( Mesh, HaloNode, CheckForHalo ) 
   END IF

   ! Go over all boundary elements and mark nodes that should be 
   ! discontinuous and nodes that should be continuous 
   DO t = 1, NoBoundElems
     
     Element => Mesh % Elements(NoBulkElems + t)
     Indexes => Element % NodeIndexes
     n = Element % Type % NumberOfNodes

     DisCont = .FALSE.
     DO bc = 1,Model % NumberOfBCs
       IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
         DisCont = ListGetLogical( Model % BCs(bc) % Values,'Discontinuous Boundary',Found )
         IF( DisCont ) EXIT
       END IF
     END DO     
     IF(.NOT. DisCont ) CYCLE
     
     DO i=1,n
       j = Indexes(i) 
       IF( CheckForHalo ) THEN
         IF( HaloNode(j) ) CYCLE
       END IF
       DisContNode(j) = .TRUE.
     END DO
     DisContElem( t ) = .TRUE.
     
     LeftElem => Element % BoundaryInfo % Left
     IF( ASSOCIATED( LeftElem ) ) THEN
       ParentUsed( LeftElem % ElementIndex ) = .TRUE.
     ELSE
       NoMissingElems = NoMissingElems + 1 
     END IF
     
     RightElem => Element % BoundaryInfo % Right
     IF( ASSOCIATED( RightElem ) ) THEN
       ParentUsed( RightElem % ElementIndex ) = .TRUE.
     ELSE
       NoMissingElems = NoMissingElems + 1
     END IF
   END DO
   
   IF( NoMissingElems > 0 ) THEN
     CALL Warn(Caller,'Missing '//TRIM(I2S(NoMissingElems))// &
     ' parent elements in partition '//TRIM(I2S(ParEnv % MyPe))) 
   END IF

   ! Calculate the number of discontinuous nodes and the number of bulk elements 
   ! associated to them. 
   NoDisContElems = COUNT( DiscontElem )
   NoDisContNodes = COUNT( DisContNode ) 
   CALL Info(Caller,'Number of discontinuous boundary elements: '&
       //TRIM(I2S(NoDisContElems)),Level=7)
   CALL Info(Caller,'Number of candicate nodes: '&
       //TRIM(I2S(NoDisContNodes)),Level=7)

   ! By default all nodes that are associated to elements immediately at the discontinuous 
   ! boundary are treated as discontinuous. However, the user may be not be greedy and release
   ! some nodes from the list that are associated also with other non-discontinuous elements.   
   ConflictElems = 0
   IF( NoDiscontNodes > 0 ) THEN
     n = NoDiscontNodes
     
     GreedyBulk = ListGetLogical( Model % Simulation,'Discontinuous Bulk Greedy',Found ) 
     IF(.NOT. Found ) GreedyBulk = .TRUE.     
     
     GreedyBC = ListGetLogical( Model % Simulation,'Discontinuous Boundary Greedy',Found ) 
     IF(.NOT. Found ) GreedyBC = .TRUE.     
     
     IF( .NOT. ( GreedyBC .AND. GreedyBulk ) ) THEN
       CALL Info(Caller,'Applying non-greedy strategies for Discontinuous mesh',Level=12)

       DO t = 1,NoBulkElems+NoBoundElems
         Element => Mesh % Elements(t)

         IF( t <= NoBulkElems ) THEN
           IF( GreedyBulk ) CYCLE
           IF( ParentUsed(t) ) CYCLE
         ELSE
           IF( GreedyBC ) CYCLE
           IF( DiscontElem(t-NoBulkElems) ) CYCLE
           !IF( Element % BoundaryInfo % Constraint == 0 ) CYCLE
           ! Check that this is not an internal BC
           IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) CYCLE
           IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right) ) CYCLE
         END IF
         Indexes => Element % NodeIndexes

         IF( ANY( DisContNode( Indexes ) ) ) THEN
           !PRINT *,'t',Element % BoundaryInfo % Constraint, t,DisContElem(t), &
           !    Indexes, DisContNode( Indexes ) 
           DisContNode( Indexes ) = .FALSE.
           ConflictElems = ConflictElems + 1
         END IF
       END DO
       NoDisContNodes = COUNT( DisContNode ) 
     END IF

     IF( ConflictElems > 0 ) THEN
       CALL Info(Caller,'Conflicting discontinuity in elements: '&
           //TRIM(I2S(ConflictElems)))
     END IF

     IF( NoDiscontNodes < n ) THEN
       CALL Info(Caller,'Number of local discontinuous nodes: '&
           //TRIM(I2S(NoDisContNodes)), Level=12)
     ELSE
       CALL Info(Caller,'All candidate nodes used',Level=12)
     END IF
     
     IF( NoDiscontNodes == 0 ) THEN
       IF( n > 0 .AND. .NOT. GreedyBulk ) THEN
         CALL Info(Caller,'You might want to try the Greedy bulk strategy',Level=3)
       END IF
     END IF
   END IF
   
   i = ParallelReduction( NoDiscontNodes ) 
   CALL Info(Caller,'Number of discontinuous nodes: '&
       //TRIM(I2S(i)),Level=7)

   IF( i == 0 ) THEN
     CALL Warn(Caller,'Nothing to create, exiting...')
     IF( CheckForHalo ) DEALLOCATE( HaloNode ) 
     DEALLOCATE( DiscontNode, DiscontElem, ParentUsed )
     RETURN
   END IF

   ! Ok, we have marked discontinuous nodes, now give them an index. 
   ! This should also create the indexes in parallel.
   DisContPerm => NULL()
   ALLOCATE( DisContPerm(NoNodes) )
   DisContPerm = 0    

   ! We could end up here on an parallel case only
   ! Then we must make the parallel numbering, so jump to the end where this is done. 
   IF( NoDisContNodes == 0 ) THEN
     IF( DoubleBC ) THEN       
       Mesh % DiscontMesh = .FALSE.
       DEALLOCATE( DisContPerm ) 
     ELSE
       Mesh % DisContMesh = .TRUE.
       Mesh % DisContPerm => DisContPerm
       Mesh % DisContNodes = 0
     END IF
     GOTO 200
   END IF
   
   ! Create a table showing nodes that are related to the moving nodes by
   ! the moving elements. 
   ALLOCATE( MovingNode( NoNodes ), StayingNode( NoNodes ) ) 
   MovingNode = .FALSE.
   StayingNode = .FALSE.

   ! For historical reasons there is both single 'body' and multiple 'bodies'
   ! that define on which side of the discontinuity the new nodes will be. 
   DiscontFlag = 'Discontinuous Target Bodies'
   TargetBodies => ListGetIntegerArray( BCList, DiscontFlag, UseTargetBodies ) 
   IF(.NOT. UseTargetBodies ) THEN
     DiscontFlag = 'Discontinuous Target Body'
     TargetBodies => ListGetIntegerArray( BCList, DiscontFlag, UseTargetBodies ) 
   END IF

   ! If either parent is consistently one of the bodies then we can create a discontinuous 
   ! boundary. Note that this currently only works currently in serial!
   IF(.NOT. UseTargetBodies ) THEN
     IF( ParEnv % PEs > 1 ) THEN
       CALL Fatal(Caller,'Please give > Discontinuous Target Bodies < on the BC!')
     END IF
     
     CALL Info(Caller,'Trying to find a dominating parent body',Level=12)

     CandA = -1
     CandB = -1
     DO t=1, NoBoundElems
       IF(.NOT. DisContElem(t) ) CYCLE
       Element => Mesh % Elements(NoBulkElems + t)

       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
         CALL Fatal(Caller,'Alternative strategy requires all parent elements!')
       END IF
       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
         CALL Fatal(Caller,'Alternative strategy requires all parent elements!')
       END IF

       LeftBody = Element % BoundaryInfo % Left % BodyId         
       RightBody = Element % BoundaryInfo % Right % BodyId

       IF( CandA == -1 ) THEN
         CandA = LeftBody 
       ELSE IF( CandA == 0 ) THEN
         CYCLE
       ELSE IF( CandA /= LeftBody .AND. CandA /= RightBody ) THEN
         CandA = 0
       END IF

       IF( CandB == -1 ) THEN
         CandB = RightBody
       ELSE IF( CandB == 0 ) THEN
         CYCLE
       ELSE IF( CandB /= LeftBody .AND. CandB /= RightBody ) THEN
         CandB = 0
       END IF
     END DO

     ! Choose the bigger one to honor the old convention
     ! This eliminates at the same time the unsuccessful case of zero.
     TargetBody(1) = MAX( CandA, CandB )

     IF( TargetBody(1) > 0 ) THEN
       CALL Info(Caller,&
           'There seems to be a consistent discontinuous body: '&
           //TRIM(I2S(TargetBody(1))),Level=8)
       UseConsistantBody = .TRUE.
       TargetBodies => TargetBody
     ELSE
       CALL Fatal(Caller,&
           'No simple rules available for determining discontinuous body')
     END IF
   END IF


   ! Assume we have only one active BC and we know the list of discontinuous 
   ! target bodies there. Hence we have all the info needed to set the 
   ! discontinuous elements also for other bulk elements. 
   ! This could be made more generic...
   NoUndecided = 0
   NoMovingElems = 0 
   NoStayingElems = 0

   DO t=1, NoBulkElems
     Element => Mesh % Elements(t)

     ! No need to treat halo elements
     !IF( CheckForHalo .AND. Element % PartIndex /= ParEnv % MyPe ) CYCLE

     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE
     Moving = ANY( TargetBodies == Element % BodyId )

     IF( Moving ) THEN
       NoMovingElems = NoMovingElems + 1 
       MovingNode(Indexes) = .TRUE.
     ELSE
       StayingNode(Indexes) = .TRUE.
       NoStayingElems = NoStayingElems + 1
     END IF
   END DO

   CALL Info(Caller,'Number of bulk elements moving: '&
       //TRIM(I2S(NoMovingElems)), Level=8)
   CALL Info(Caller,'Number of bulk elements staying: '&
       //TRIM(I2S(NoStayingElems)), Level=8)

   ! Set discontinuous nodes only if there is a real moving node associated with it
   ! Otherwise we would create a zero to the permutation vector. 
   ! If there is just a staying node then no need to create discontinuity at this node.
   DiscontNode = DiscontNode .AND. MovingNode 

   ! Create permutation numbering for the discontinuous nodes   
   ! Doubling will be done only for nodes that have both parents
   j = 0
   DO i=1,NoNodes
     IF( DisContNode(i) ) THEN
       j = j + 1
       DisContPerm(i) = j
     END IF
   END DO
   IF( j < NoDiscontNodes ) THEN
     PRINT *,'Some discontinuous nodes only needed on the other side:',&
         ParEnv % MyPe, NoDiscontNodes-j
     NoDiscontNodes = j 
   END IF


   ! Now set the new indexes for bulk elements
   ! In parallel skip the halo elements
   DO t=1, NoBulkElems
     Element => Mesh % Elements(t)

     ! No need to treat halo elements
     !IF( CheckForHalo .AND. Element % PartIndex /= ParEnv % MyPe ) CYCLE
     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE
     Moving = ANY( TargetBodies == Element % BodyId )

     IF( Moving ) THEN
       DO i=1, SIZE(Indexes) 
         j = DisContPerm(Indexes(i))
         IF( j > 0 ) Indexes(i) = NoNodes + j
       END DO
     END IF
   END DO

    
   ! Now set also the unset boundary elements by following the ownership of the parent elements
   ! or the majority opinion if this is conflicting.
   DO t=1, NoBoundElems

     Element => Mesh % Elements(NoBulkElems + t)

     ! If the element has no constraint then there is no need to treat it
     IF( Element % BoundaryInfo % Constraint == 0 ) CYCLE

     IF( DisContElem(t) ) THEN
       LeftElem => Element % BoundaryInfo % Left
       RightElem => Element % BoundaryInfo % Right

       IF( ASSOCIATED( LeftElem ) ) THEN
         Moving = ANY( TargetBodies == LeftElem % BodyId ) 
       ELSE
         Moving = .NOT. ANY( TargetBodies == RightElem % BodyId )
       END IF
       IF( Moving ) THEN
         Element % BoundaryInfo % Left => RightElem
         Element % BoundaryInfo % Right => LeftElem 
       END IF
       CYCLE
     END IF


     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE

     ElemFamily = Element % TYPE % ElementCode / 100 
     LeftElem => Element % BoundaryInfo % Left
     RightElem => Element % BoundaryInfo % Right

     ! The boundary element follows the parent element if it is clear what to do
     Set = .TRUE.
     IF( ASSOCIATED( LeftElem ) .AND. ASSOCIATED( RightElem ) ) THEN
       Moving = ANY( TargetBodies == LeftElem % BodyId )
       Moving2 = ANY( TargetBodies == RightElem % BodyId ) 
       IF( Moving .NEQV. Moving2) THEN
         CALL Warn(Caller,'Conflicting moving information')
         !PRINT *,'Moving:',t,Element % BoundaryInfo % Constraint, &
         !    Moving,Moving2,LeftElem % BodyId, RightElem % BodyId
         Set = .FALSE.
       ELSE
         IF( Moving ) THEN
           Element % BoundaryInfo % Left => RightElem
           Element % BoundaryInfo % Right => LeftElem 
         END IF
       END IF
     ELSE IF( ASSOCIATED( LeftElem ) ) THEN
       Moving = ANY( LeftElem % NodeIndexes > NoNodes ) 
     ELSE IF( ASSOCIATED( RightElem ) ) THEN
       Moving = ANY( RightElem % NodeIndexes > NoNodes )
     ELSE
       CALL Fatal(Caller,'Boundary BC has no parants!')
     END IF

     ! Otherwise we follow the majority rule
     IF( .NOT. Set ) THEN
       NoMoving = COUNT( MovingNode(Indexes) ) 
       NoStaying = COUNT( StayingNode(Indexes) ) 

       IF( NoStaying /= NoMoving ) THEN
         Moving = ( NoMoving > NoStaying )
         Set = .TRUE.
       END IF
     END IF

     ! Ok, finally set whether boundary element is moving or staying
     IF( Set ) THEN
       IF( Moving ) THEN
         NoMovingElems = NoMovingElems + 1 
         DO i=1, SIZE(Indexes) 
           j = DisContPerm(Indexes(i))
           IF( j > 0 ) Indexes(i) = NoNodes + j
         END DO
       ELSE
         NoStayingElems = NoStayingElems + 1
       END IF
     ELSE
       NoUndecided = NoUndecided + 1
     END IF
   END DO

   CALL Info(Caller,'Number of related elements moving: '&
       //TRIM(I2S(NoMovingElems)), Level=8 )
   CALL Info(Caller,'Number of related elements staying: '&
       //TRIM(I2S(NoStayingElems)), Level=8 )
   IF( NoUndecided == 0 ) THEN
     CALL Info(Caller,'All elements marked either moving or staying')
   ELSE
     CALL Info(Caller,'Number of related undecided elements: '//TRIM(I2S(NoUndecided)) )
     CALL Warn(Caller,'Could not decide what to do with some boundary elements!')
   END IF


   m = COUNT( DiscontNode .AND. .NOT. MovingNode )
   IF( m > 0 ) THEN
     PRINT *,'Number of discont nodes not moving: ',ParEnv % MyPe, m
   END IF

   m = COUNT( DiscontNode .AND. .NOT. StayingNode )
   IF( m > 0 ) THEN
     PRINT *,'Number of discont nodes not staying: ',ParEnv % MyPe, m
     DO i=1,SIZE(DisContNode)
       IF( DiscontNode(i) .AND. .NOT. StayingNode(i) ) THEN
         IF( ParEnv % PEs == 1 ) THEN
           PRINT *,'Node:',ParEnv % MyPe,i
         ELSE
           PRINT *,'Node:',ParEnv % MyPe,i,Mesh % ParallelInfo % GlobalDofs(i), &
               Mesh % ParallelInfo % NeighbourList(i) % Neighbours
         END IF
         PRINT *,'Coord:',ParEnv % MyPe, Mesh % Nodes % x(i), Mesh % Nodes % y(i)
       END IF
     END DO
   END IF

   !DEALLOCATE( MovingNode, StayingNode )

   ! Now add the new nodes also to the nodes structure
   ! and give the new nodes the same coordinates as the ones
   ! that they were derived from. 
   Mesh % NumberOfNodes = NoNodes + NoDisContNodes   
   CALL EnlargeCoordinates( Mesh ) 

   CALL Info(Caller,'Setting new coordinate positions',Level=12)
   DO i=1, NoNodes
     j = DisContPerm(i)
     IF( j > 0 ) THEN
       k = NoNodes + j
       Mesh % Nodes % x(k) = Mesh % Nodes % x(i)
       Mesh % Nodes % y(k) = Mesh % Nodes % y(i)
       Mesh % Nodes % z(k) = Mesh % Nodes % z(i)
     END IF
   END DO


   ! If the discontinuous boundary is duplicated then no information of it 
   ! is saved. The periodic and mortar conditions now need to perform
   ! searches. On the other hand the meshes may now freely move.,
   IF( DoubleBC ) THEN
     CALL Info(Caller,'Creating secondary boundary for Discontinuous gap',Level=10)

     CALL EnlargeBoundaryElements( Mesh, NoDiscontElems ) 

     NoDisContElems = 0
     DO t=1, NoBoundElems

       ! Is this a boundary to be doubled?
       IF(.NOT. DisContElem(t) ) CYCLE

       Element => Mesh % Elements(NoBulkElems + t)
       IF(.NOT. ASSOCIATED(Element) ) THEN
         CALL Fatal(Caller,'Element '//TRIM(I2S(NoBulkElems+t))//' not associated!')
       END IF
       Indexes => Element % NodeIndexes

       DisContTarget = 0
       Found = .FALSE.
       DO bc = 1,Model % NumberOfBCs
         IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Discontinuous BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Mortar BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Periodic BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Contact BC',Found )
           IF( Found ) EXIT
         END IF
       END DO
       IF( .NOT. Found .OR. DisContTarget == 0 ) THEN
         CALL Fatal(Caller,'Nonzero target boundary must be given for all, if any bc!')
       END IF

       RightElem => Element % BoundaryInfo % Right
       LeftElem => Element % BoundaryInfo % Left 

       NoDisContElems = NoDisContElems + 1              
       j = NoBulkElems + NoBoundElems + NoDisContElems 

       OtherElem => Mesh % Elements( j )
       IF(.NOT. ASSOCIATED(OtherElem) ) THEN
         CALL Fatal(Caller,'Other elem '//TRIM(I2S(j))//' not associated!')
       END IF

       OtherElem = Element 
       OtherElem % TYPE => Element % TYPE

       NULLIFY( OtherElem % BoundaryInfo ) 
       ALLOCATE( OtherElem % BoundaryInfo ) 
       OtherElem % BoundaryInfo % Left => Element % BoundaryInfo % Right

       ! Now both boundary elements are just one sided. Remove the associated to the other side. 
       NULLIFY( Element % BoundaryInfo % Right ) 
       NULLIFY( OtherElem % BoundaryInfo % Right )

       NULLIFY( OtherElem % NodeIndexes )
       n = SIZE( Element % NodeIndexes ) 
       ALLOCATE( OtherElem % NodeIndexes( n ) ) 

       ! Ok, we found the element to manipulate the indexes. 
       ! The new index is numbered on top of the old indexes. 
       DO i=1,n
         j = Element % NodeIndexes(i) 
         IF( DisContPerm(j) > 0 ) THEN
           OtherElem % NodeIndexes(i) = NoNodes + DisContPerm(j)
         ELSE 
           OtherElem % NodeIndexes(i) = j
         END IF
       END DO

       OtherElem % BoundaryInfo % Constraint = DisContTarget
     END DO

     CALL Info(Caller,'Number of original bulk elements: '&
         //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=10)
     CALL Info(Caller,'Number of original boundary elements: '&
         //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=10)
     CALL Info(Caller,'Number of additional boundary elements: '&
         //TRIM(I2S(NoDisContElems)),Level=10)

     Mesh % DiscontMesh = .FALSE.
   ELSE
     Mesh % DisContMesh = .TRUE.
     Mesh % DisContPerm => DisContPerm
     Mesh % DisContNodes = NoDisContNodes 
   END IF

200 CONTINUE


   CALL EnlargeParallelInfo(Mesh, DiscontPerm )
   IF( ParEnv % PEs > 1 ) THEN
     m = COUNT( Mesh % ParallelInfo % GlobalDofs == 0) 
     IF( m > 0 ) CALL Warn(Caller,'There are nodes with zero global dof index: '//TRIM(I2S(m)))
   END IF

   IF( DoubleBC .AND. NoDiscontNodes > 0 ) DEALLOCATE( DisContPerm )


   DEALLOCATE( DisContNode, DiscontElem )   
  
 END SUBROUTINE CreateDiscontMesh


!> Reallocate coordinate arrays for iso-parametric p-elements,
!> or if the size of nodes has been increased due to discontinuity. 
!> This does not seem to be necessary for other types of 
!> elements (face, edge, etc.)
! -----------------------------------------------------------    
 SUBROUTINE EnlargeCoordinates(Mesh)

   TYPE(Mesh_t) :: Mesh
   INTEGER :: n0, n
   REAL(KIND=dp), POINTER :: TmpCoord(:)

   INTEGER :: i
   LOGICAL :: pelementsPresent

   n = Mesh % NumberOfNodes + &
       Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
       Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
       Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements
   n0 = SIZE( Mesh % Nodes % x )

   pelementsPresent = .FALSE.
   DO i=1,Mesh % NumberOfBulkElements
     IF(isPelement(Mesh % Elements(i))) THEN
       pelementsPresent = .TRUE.; EXIT
     END IF
   END DO

   IF ( Mesh % NumberOfNodes > n0 .OR. n > n0 .AND. pelementsPresent ) THEN
     CALL Info('EnlargeCoordinates','Increasing number of nodes from '&
         //TRIM(I2S(n0))//' to '//TRIM(I2S(n)),Level=8)

     TmpCoord => Mesh % Nodes % x
     ALLOCATE( Mesh % Nodes % x(n) )
     Mesh % Nodes % x(1:n0) = TmpCoord
     Mesh % Nodes % x(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )

     TmpCoord => Mesh % Nodes % y
     ALLOCATE( Mesh % Nodes % y(n) )
     Mesh % Nodes % y(1:n0) = TmpCoord
     Mesh % Nodes % y(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )

     TmpCoord => Mesh % Nodes % z
     ALLOCATE( Mesh % Nodes % z(n) )
     Mesh % Nodes % z(1:n0) = TmpCoord
     Mesh % Nodes % z(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )
   END IF

 END SUBROUTINE EnlargeCoordinates


 
 SUBROUTINE EnlargeBoundaryElements(Mesh, DoubleElements )

   TYPE(Mesh_t) :: Mesh
   INTEGER :: DoubleElements
   INTEGER :: n,n0,i,j
   REAL(KIND=dp), POINTER :: TmpCoord(:)
   TYPE(Element_t), POINTER :: NewElements(:),OldElements(:), Element

   IF( DoubleElements == 0 ) RETURN

   n0 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
   n = n0 + DoubleElements

   CALL Info('EnlargeBoundaryElements','Increasing number of elements from '&
       //TRIM(I2S(n0))//' to '//TRIM(I2S(n)),Level=8)

   OldElements => Mesh % Elements
   CALL AllocateVector( Mesh % Elements, n, 'EnlargeBoundaryElements' )
   DO i=1,n0
     Mesh % Elements(i) = OldElements(i)
     IF(ASSOCIATED(OldElements(i) % BoundaryInfo)) THEN
       IF (ASSOCIATED(OldElements(i) % BoundaryInfo % Left)) &
           Mesh % Elements(i) % BoundaryInfo % Left => &
           Mesh % Elements(OldElements(i) % BoundaryInfo % Left % ElementIndex)
       
       IF (ASSOCIATED(OldElements(i) % BoundaryInfo % Right)) &
           Mesh % Elements(i) % BoundaryInfo % Right => &
           Mesh % Elements(OldElements(i) % BoundaryInfo % Right % ElementIndex)
     END IF
   END DO

   DO i=n0+1,n
     Element => Mesh % Elements(i)

     Element % DGDOFs = 0
     Element % BodyId = 0
     Element % TYPE => NULL()
     Element % BoundaryInfo => NULL()
     Element % PDefs => NULL()
     Element % DGIndexes => NULL()
     Element % EdgeIndexes => NULL()
     Element % FaceIndexes => NULL()
     Element % BubbleIndexes => NULL()
   END DO

   DEALLOCATE( OldElements ) 
   Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements + DoubleElements

 END SUBROUTINE EnlargeBoundaryElements


 SUBROUTINE EnlargeParallelInfo( Mesh, DiscontPerm )

   TYPE(Mesh_t) :: Mesh
   INTEGER, POINTER :: DiscontPerm(:)

   INTEGER :: nmax,n0,n1,i,j,istat, goffset
   INTEGER, POINTER :: TmpGlobalDofs(:) 
   INTEGER, ALLOCATABLE :: Perm(:)
   LOGICAL, POINTER :: Intf(:)
   TYPE(NeighbourList_t), POINTER :: Nlist(:)

   IF ( ParEnv % PEs <= 1 ) RETURN

   ! As index offset use the number of nodes in the whole mesh
   goffset = ParallelReduction( MAXVAL(Mesh % ParallelInfo % GlobalDofs),2 )

   n0 = SIZE( Mesh % ParallelInfo % GlobalDofs )
   n1 = Mesh % NumberOfNodes 
   IF( n0 >= n1 ) THEN
     CALL Info('EnlargeParallelInfo','No need to grow: '&
         //TRIM(I2S(n0))//' vs. '//TRIM(I2S(n1)),Level=10)
     RETURN
   END IF
   
   CALL Info('EnlargeParallelInfo','Increasing global numbering size from '&
         //TRIM(I2S(n0))//' to '//TRIM(I2S(n1)),Level=8)

   ! Create permutation table for the added nodes
   ALLOCATE(Perm(n1)); Perm  = 0
   DO i=1,n0
     IF ( DiscontPerm(i) > 0 ) THEN
       Perm(DiscontPerm(i)+n0) = i
     END IF
   END DO

   ! Create the enlarged set of global nodes indexes
   ALLOCATE( TmpGlobalDofs(n1), STAT=istat )
   IF (istat /= 0) CALL Fatal('EnlargeParallelInfo', 'Unable to allocate TmpGlobalDofs array.')
   TmpGlobalDofs = 0
   DO i=1,n0
     TmpGlobalDofs(i) = Mesh % ParallelInfo % GlobalDofs(i)
   END DO
   DO i=n0+1,n1
     j = Perm(i)
     IF(j > 0) THEN
       TmpGlobalDofs(i) = TmpGlobalDOfs(j) + goffset
     END IF
   END DO
   DEALLOCATE(Mesh % ParallelInfo % GlobalDofs)
   Mesh % ParallelInfo % GlobalDOfs => TmpGlobalDofs

   ! Create the enlarged list of neighbours
   ALLOCATE(Nlist(n1))
   DO i=1,n0
     IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
       Nlist(i) % Neighbours => &
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       Mesh % ParallelInfo % NeighbourList(i) % Neighbours => NULL()
     ELSE 
       Nlist(i) % Neighbours => NULL()
     END IF
   END DO

   DO i=n0+1,n1
     j = Perm(i)
     IF ( j > 0 ) THEN
       IF( ASSOCIATED( Nlist(j) % Neighbours ) ) THEN
         ALLOCATE( Nlist(i) % Neighbours(SIZE(Nlist(j) % Neighbours) ) )
         Nlist(i) % Neighbours = Nlist(j) % Neighbours
       ELSE
         Nlist(i) % Neighbours => NULL()
       END IF
     END IF
   END DO
   DEALLOCATE(Mesh % ParallelInfo % NeighbourList)
   Mesh % ParallelInfo % NeighbourList => Nlist


   ! Create logical table showing the interface nodes
   ALLOCATE( Intf(n1) )
   Intf = .FALSE.
   Intf(1:n0) = Mesh % ParallelInfo % NodeInterface(1:n0)
   DO i=n0+1,n1
     j = Perm(i)
     IF(j > 0 ) THEN
       Intf(i) = Intf(j) 
     END IF
   END DO
   DEALLOCATE( Mesh % ParallelInfo % NodeInterface )
   Mesh % ParallelInfo % NodeInterface => Intf


 END SUBROUTINE EnlargeParallelInfo




 !> Fortran reader for Elmer ascii mesh file format.
 !> This is a Fortran replacement for the old C++ eio library. 
 !------------------------------------------------------------------------
 SUBROUTINE ElmerAsciiMesh(Step, PMesh, MeshNamePar, ThisPe, NumPEs, IsParallel )

   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel

   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: PrevStep=0, iostat
   INTEGER, PARAMETER :: FileUnit = 10
   CHARACTER(MAX_NAME_LEN) :: BaseName, FileName
   INTEGER :: i,j,k,n,BaseNameLen, SharedNodes = 0, mype = 0, numprocs = 0
   INTEGER, POINTER :: NodeTags(:), ElementTags(:), LocalPerm(:)
   INTEGER :: MinNodeTag = 0, MaxNodeTag = 0, istat
   LOGICAL :: ElementPermutation=.FALSE., NodePermutation=.FALSE., Parallel, &
       PseudoParallel, Found


   SAVE PrevStep, BaseName, BaseNameLen, Mesh, mype, Parallel, &
       NodeTags, ElementTags, LocalPerm, PseudoParallel

   CALL Info('ElmerAsciiMesh','Performing step: '//TRIM(I2S(Step)),Level=8)

   IF( Step - PrevStep /= 1 ) THEN
     CALL Fatal('ElmerAsciiMesh','The routine should be called in sequence: '// &
         TRIM(I2S(PrevStep))//' : '//TRIM(I2S(Step)) )
   END IF
   PrevStep = Step
   IF( PrevStep == 6 ) PrevStep = 0 

   IF( Step == 1 ) THEN
     IF(.NOT. PRESENT( MeshNamePar ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give MeshNamePar!')
     END IF
     BaseName = TRIM( MeshNamePar ) 
     IF(.NOT. PRESENT( PMesh ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give PMesh!')
     END IF
     Mesh => PMesh
     IF(.NOT. PRESENT( ThisPe ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give ThisPe!')
     END IF
     mype = ThisPe 
     IF(.NOT. PRESENT( NumPEs) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give NumPEs!')
     END IF
     numprocs = NumPEs
     IF(.NOT. PRESENT( IsParallel ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give IsParallel!')
     END IF
     Parallel = IsParallel

     PseudoParallel = .FALSE.
     IF(.NOT. Parallel ) THEN
       IF( ParEnv % PEs > 1 ) THEN
         PseudoParallel = ListGetLogical(CurrentModel % Simulation,'Enforce Parallel',Found ) 
         IF(.NOT. Found ) PseudoParallel = ListGetLogicalAnySolver(CurrentModel,'Enforce Parallel')
       END IF
     END IF
     
     i = LEN_TRIM(MeshNamePar)
     DO WHILE(MeshNamePar(i:i) == CHAR(0))
       i=i-1
     END DO
     BaseNameLen = i
     CALL Info('ElmerAsciiMesh','Base mesh name: '//TRIM(MeshNamePar(1:BaseNameLen)))
   END IF
   

   SELECT CASE( Step ) 

   CASE(1)       
     CALL ReadHeaderFile()

   CASE(2)
     CALL ReadNodesFile()

   CASE(3)
     CALL ReadElementsFile()

   CASE(4)
     CALL ReadBoundaryFile()
     CALL PermuteNodeNumbering()

   CASE(5)
     IF( PseudoParallel ) THEN
       CALL InitPseudoParallel()
     ELSE
       CALL InitParallelInfo()
       CALL ReadSharedFile()
     END IF
       
   CASE(6)
     IF( ASSOCIATED( LocalPerm) ) DEALLOCATE( LocalPerm ) 
     IF( ASSOCIATED( ElementTags) ) DEALLOCATE( ElementTags )

   END SELECT


 CONTAINS


   FUNCTION read_ints(s,j,halo) RESULT(n)
     INTEGER :: j(:)
     CHARACTER(LEN=*) :: s
     LOGICAL :: halo
     
     INTEGER :: i,k,l,m,n,ic
     INTEGER, PARAMETER :: ic0 = ICHAR('0'), ic9 = ICHAR('9'), icm = ICHAR('-'), &
         icd = ICHAR('/'), ics = ICHAR(' ')
     
     k = LEN_TRIM(s)
     l = 1
     n = 0
     halo = .FALSE.
     DO WHILE(l<=k.AND.n<SIZE(j))
       DO WHILE(l<=k)
         ic = ICHAR(s(l:l))
         IF( ic == ics ) THEN
           CONTINUE
         ELSE IF( ic == icd ) THEN
           halo = .TRUE.
         ELSE
           EXIT
         END IF
         l=l+1
       END DO
       IF(l>k) EXIT
       IF(.NOT.(ic==icm .OR. ic>=ic0 .AND. ic<=ic9)) EXIT
       
       m = l+1
       DO WHILE(m<=k)
         ic = ICHAR(s(m:m))
         IF(ic<ic0 .OR. ic>ic9) EXIT
         m=m+1
       END DO
       
       n = n + 1
       j(n) = s2i(s(l:m-1),m-l)
       l = m
     END DO
   END FUNCTION read_ints
   

   !---------------------------------------------------
   ! Read header file and allocate some mesh structures
   !---------------------------------------------------
   SUBROUTINE ReadHeaderFile()

     INTEGER :: TypeCount
     INTEGER :: Types(64),CountByType(64)

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(numprocs))//&
           '/part.'//TRIM(I2S(mype+1))//'.header'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.header'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadHeaderFile','Reading header info from file: '//TRIM(FileName),Level=10)
     END IF

     READ(FileUnit,*,IOSTAT=iostat) Mesh % NumberOfNodes, &
         Mesh % NumberOfBulkElements,&
         Mesh % NumberOfBoundaryElements
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not read header 1st line in file: '//TRIM(FileName))
     END IF

     Types = 0
     CountByType = 0
     READ(FileUnit,*,IOSTAT=iostat) TypeCount
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not read the type count in file: '//TRIM(FileName))
     END IF
     DO i=1,TypeCount
       READ(FileUnit,*,IOSTAT=iostat) Types(i),CountByType(i)
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadHeaderFile','Could not read type count '&
             //TRIM(I2S(i))//'in file: '//TRIM(FileName))
       END IF
     END DO

     IF( Parallel ) THEN
       READ(FileUnit,*,IOSTAT=iostat) SharedNodes
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadHeaderFile','Could not read shared nodes in file: '//TRIM(FileName))
       END IF
     ELSE
       SharedNodes = 0
     END IF

     Mesh % MaxElementNodes = 0
     DO i=1,TypeCount
       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes, MODULO( Types(i), 100) )
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadHeaderFile


   !-----------------------------------------------------------------------
   ! Read nodes file and create nodal permutation if needed
   !-----------------------------------------------------------------------
   SUBROUTINE ReadNodesFile()

     REAL(KIND=dp) :: Coords(3)
     INTEGER :: NodeTag

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(numprocs))//&
           '/part.'//TRIM(I2S(mype+1))//'.nodes'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.nodes'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadNodesFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadNodesFile','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ALLOCATE( NodeTags(Mesh % NumberOfNodes ) ) 
     NodeTags = 0

     NodePermutation = .FALSE.
     DO j = 1, Mesh % NumberOfNodes
       READ(FileUnit,*,IOSTAT=iostat) NodeTag, k, Coords
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadNodesFile','Problem load node '//TRIM(I2S(j))//' in file: '//TRIM(Filename))
       END IF

       IF( NodeTags(j) /= j ) NodePermutation = .TRUE.
 
       NodeTags(j) = NodeTag
       Mesh % Nodes % x(j) = Coords(1)
       Mesh % Nodes % y(j) = Coords(2)
       Mesh % Nodes % z(j) = Coords(3)
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadNodesFile


   !------------------------------------------------------------------------------
   ! Read elements file and create elemental permutation if needed 
   !------------------------------------------------------------------------------
   SUBROUTINE ReadElementsFile()
     TYPE(Element_t), POINTER :: Element
     INTEGER :: ElemType, Tag, Body, ElemNo, Ivals(64),nread, ioffset, partn
     CHARACTER(256) :: str
     LOGICAL :: halo


     CALL AllocateVector( ElementTags, Mesh % NumberOfBulkElements+1, 'ReadElementsFile')   
     ElementTags = 0
     ElementPermutation = .FALSE.

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)// &
          '/partitioning.'//TRIM(I2S(numprocs))//&
             '/part.'//TRIM(I2S(mype+1))//'.elements'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.elements'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', iostat=IOSTAT )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadElementsFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadElementsFile','Reading bulk elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=1,Mesh % NumberOfBulkElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadElementsFile','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadElementsFile','Could not read start of element entry: '//TRIM(I2S(j)))
       END IF

       nread = read_ints(str,ivals,halo)

       tag = ivals(1)

       IF( halo ) THEN
         ioffset = 1
         partn = ivals(2) 
       ELSE
         ioffset = 0
         partn = 0 
       END IF
       body = ivals(ioffset+2)
       ElemType = ivals(ioffset+3)

       ElementTags(j) = tag
       IF( j /= tag ) ElementPermutation = .TRUE.             
       Element % ElementIndex = j
       Element % BodyId = body

       IF( partn > 0 ) THEN
         Element % PartIndex = partn-1
       ELSE
         Element % PartIndex = mype
       END IF

       Element % TYPE => GetElementType(ElemType)

       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadElementsFile','Element of type '&
             //TRIM(I2S(ElemType))//' could not be associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       IF( nread < n + ioffset + 3 ) THEN
         CALL Fatal('ReadElementsFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF

       CALL AllocateVector( Element % NodeIndexes, n )

       Element % NodeIndexes(1:n) = IVals(4+ioffset:nread)
     END DO
     CLOSE( FileUnit ) 

   END SUBROUTINE ReadElementsFile
   !------------------------------------------------------------------------------


   !------------------------------------------------------------------------------
   ! Read boundary elements file and remap the parents if needed.  
   !------------------------------------------------------------------------------
   SUBROUTINE ReadBoundaryFile()
     INTEGER, POINTER :: LocalEPerm(:)
     INTEGER :: MinEIndex, MaxEIndex, ElemNodes, i
     INTEGER :: Left, Right, bndry, tag, ElemType, IVals(64), nread, ioffset, partn
     TYPE(Element_t), POINTER :: Element
     CHARACTER(256) :: str
     LOGICAL :: halo

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(numprocs))//&
           '/part.'//TRIM(I2S(mype+1))//'.boundary'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.boundary'
     END IF

     ! Create permutation for the elements. This is needed when the element 
     ! parents are mapped to the new order. This is needed for mapping of the 
     ! parents. Otherwise the element numbering is arbitrary. 
     !------------------------------------------------------------------------------
     IF( ElementPermutation ) THEN
       MinEIndex = MINVAL( ElementTags(1:Mesh % NumberOfBulkElements) )
       MaxEIndex = MAXVAL( ElementTags(1:Mesh % NumberOfBulkElements) )

       LocalEPerm => NULL()
       CALL AllocateVector( LocalEPerm, MaxEIndex - MinEIndex + 1, 'ReadBoundaryFile' )
       LocalEPerm = 0
       DO i=1,Mesh % NumberOfBulkElements
         LocalEPerm( ElementTags(i) - MinEIndex + 1 ) = i
       END DO
     ELSE
       MinEIndex = 1 
       MaxEIndex = Mesh % NumberOfBulkElements
     END IF


     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', iostat=IOSTAT )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadBoundaryFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadBoundaryFile','Reading boundary elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadBoundaryFile','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadBoundaryFile','Could not read boundary element entry: '//TRIM(I2S(j)))
       END IF
       nread = read_ints(str,ivals,halo)
       
       tag = ivals(1)

       IF( halo ) THEN
         partn = ivals(2)
         ioffset = 1
       ELSE
         partn = 0
         ioffset = 0
       END IF

       bndry = ivals(ioffset+2)
       left = ivals(ioffset+3)
       right = ivals(ioffset+4)
       ElemType = ivals(ioffset+5)
       
       Element % ElementIndex = j
       Element % TYPE => GetElementType(ElemType)
       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadBoundaryFile','Element of type '//TRIM(I2S(ElemType))//'could not be associated!')
       END IF

       ElemNodes = Element % TYPE % NumberOfNodes
       Mesh % MaxElementNodes = MAX( Mesh % MaxElementNodes, ElemNodes )

       IF( partn == 0 ) THEN
         Element % PartIndex = mype
       ELSE
         Element % PartIndex = partn-1
       END IF

       CALL AllocateBoundaryInfo( Element ) 

       Element % BoundaryInfo % Constraint = bndry
       Element % BoundaryInfo % Left => NULL()
       Element % BoundaryInfo % Right => NULL()

       IF ( Left >= MinEIndex .AND. Left <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Left  = LocalEPerm(Left - MinEIndex + 1)
         END IF
       ELSE IF ( Left > 0 ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag, Left
         CALL Error( 'ReadBoundaryFile', Message )
         Left = 0
       END IF

       IF ( Right >= MinEIndex .AND. Right <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Right = LocalEPerm(Right - MinEIndex + 1)
         END IF
       ELSE IF ( Right > 0 ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag,Right
         CALL Error( 'ReadBoundaryFile', Message )
         Right = 0
       END IF

       IF ( Left >= 1 ) THEN
         Element % BoundaryInfo % Left => Mesh % Elements(left)
       END IF

       IF ( Right >= 1 ) THEN
         Element % BoundaryInfo % Right => Mesh % Elements(right)
       END IF

       n = Element % TYPE % NumberOfNodes
       CALL AllocateVector( Element % NodeIndexes, n )

       IF( nread < 5 + n + ioffset ) THEN
         CALL Fatal('ReadBoundaryFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF
       Element % NodeIndexes(1:n) = Ivals(6+ioffset:nread)
     END DO
     CLOSE( FileUnit )


     IF( ElementPermutation ) THEN
       DEALLOCATE( LocalEPerm ) 
     END IF

   END SUBROUTINE ReadBoundaryFile
   !------------------------------------------------------------------------------



   ! Make a permutation for the bulk and boundary element topology if 
   ! the nodes are permuted. This is always the case in parallel.
   ! The initial numbering is needed only when the nodes are loaded and 
   ! hence this is a local subroutine. 
   !----------------------------------------------------------------------
   SUBROUTINE PermuteNodeNumbering()

     TYPE(Element_t), POINTER :: Element

     IF( NodePermutation ) THEN
       CALL Info('PermuteNodeNumbering','Performing node mapping',Level=6)

       MinNodeTag = MINVAL( NodeTags )
       MaxNodeTag = MAXVAL( NodeTags )

       CALL AllocateVector( LocalPerm, MaxNodeTag-MinNodeTag+1, 'PermuteNodeNumbering' )
       LocalPerm = 0
       DO i=1,Mesh % NumberOfNodes
         LocalPerm(NodeTags(i) - MinNodeTag + 1) = i
       END DO

       DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements       
         Element => Mesh % Elements(i)
         n = Element % TYPE % NumberOfNodes

         DO j=1,n
           k = Element % NodeIndexes(j) 
           Element % NodeIndexes(j) = LocalPerm(k - MinNodeTag + 1)
         END DO
       END DO
     ELSE
       CALL Info('PermuteNodeNumbering','Node mapping is continuous',Level=8)
     END IF

     ! Set the for now, if the case is truly parallel we'll have to revisit these
     ! when reading the parallel information. 
     Mesh % ParallelInfo % NumberOfIfDOFs = 0
     Mesh % ParallelInfo % GlobalDOFs => NodeTags

   END SUBROUTINE PermuteNodeNumbering


   ! Initialize some parallel structures once the non-nodal 
   ! element types are known. 
   ! Currently this is here mainly because the 
   ! Elemental and Nodal tags are local
   !-------------------------------------------------------
   SUBROUTINE InitParallelInfo()

     INTEGER, POINTER :: TmpGlobalDofs(:)

     ! These two have already been set, and if the case is serial
     ! case they can be as is.
     !Mesh % ParallelInfo % NumberOfIfDOFs = 0
     !Mesh % ParallelInfo % GlobalDOFs => NodeTags


     ! This also for serial runs ...
     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i)
     END DO

     IF(.NOT. Parallel ) RETURN

     n = Mesh % NumberOfNodes + &
         Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
         Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
         Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements

     ALLOCATE( TmpGlobalDOFs(n) )
     TmpGlobalDOFs = 0
     TmpGlobalDOFs(1:Mesh % NumberOfNodes) = &
         Mesh % ParallelInfo % GlobalDOFs(1:Mesh % NumberOfNodes)
     DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs ) 
     Mesh % ParallelInfo % GlobalDofs => TmpGlobalDofs

     ALLOCATE(Mesh % ParallelInfo % NeighbourList(n), STAT=istat)
     IF (istat /= 0) CALL Fatal('InitParallelInfo', 'Unable to allocate NeighbourList array.')

     DO i=1,n
       NULLIFY( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % NodeInterface, n, 'InitParallelInfo')
     Mesh % ParallelInfo % NodeInterface = .FALSE.       

   END SUBROUTINE InitParallelInfo


   ! Read the file that shows the shared nodes.
   !------------------------------------------------------------------------
   SUBROUTINE ReadSharedFile()

     INTEGER :: Ivals(64)
     INTEGER :: npart, tag, nread
     CHARACTER(256) :: str
     LOGICAL :: halo

     IF(.NOT. Parallel) RETURN

     FileName = BaseName(1:BaseNameLen)//&
       '/partitioning.'//TRIM(I2S(numprocs))//&
         '/part.'//TRIM(I2S(mype+1))//'.shared'

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadSharedFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadSharedFile','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ! This loop could be made more effective, for example
     ! by reading tags and nparts to a temporal vector
     ! The operation using the str takes much more time.
     !-----------------------------------------------------
     DO i=1,SharedNodes          
       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadSharedFile','Could not read shared nodes entry: '//TRIM(I2S(i)))
       END IF
       nread = read_ints(str,ivals,halo)

       tag = ivals(1)
       npart = ivals(2)       

       k = LocalPerm( tag-MinNodeTag+1 )
       Mesh % ParallelInfo % NodeInterface(k) = .TRUE.
       CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(k) % Neighbours,npart)

       IF( nread < 2 + npart ) THEN
         CALL Fatal('ReadSharedFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF
       
       Mesh % ParallelInfo % NeighbourList(k) % Neighbours = ivals(3:nread) - 1

       ! this partition does not own the node
       IF ( ivals(3)-1 /= mype ) THEN
         Mesh % ParallelInfo % NumberOfIfDOFs = &
             Mesh % ParallelInfo % NumberOfIfDOFs + 1
       END IF
     END DO

     CLOSE( FileUnit )

   END SUBROUTINE ReadSharedFile


   ! Initialize parallel info for pseudo parallel meshes
   !-------------------------------------------------------
   SUBROUTINE InitPseudoParallel()

     INTEGER, POINTER :: TmpGlobalDofs(:)

     ! This also for serial runs ...
     n = ParEnv % MyPe * Mesh % NumberOfBulkElements

     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i) + n
     END DO

     n = Mesh % NumberOfNodes + &
         Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
         Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
         Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements

     ALLOCATE( TmpGlobalDOFs(n) )
     TmpGlobalDOFs = 0
     TmpGlobalDOFs(1:Mesh % NumberOfNodes) = &
         Mesh % ParallelInfo % GlobalDOFs(1:Mesh % NumberOfNodes) + n
     DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs ) 
     Mesh % ParallelInfo % GlobalDofs => TmpGlobalDofs
     
     ALLOCATE(Mesh % ParallelInfo % NeighbourList(n), STAT=istat)
     IF (istat /= 0) CALL Fatal('InitParallelInfo', 'Unable to allocate NeighbourList array.')
     
     DO i=1,n
       ALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) )
       Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % MyPe
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % NodeInterface, n, 'InitParallelInfo')
     Mesh % ParallelInfo % NodeInterface = .FALSE.       

   END SUBROUTINE InitPseudoParallel

   
 END SUBROUTINE ElmerAsciiMesh



 !> An interface over potential mesh loading strategies. 
 !----------------------------------------------------------------- 
 SUBROUTINE LoadMeshStep( Step, PMesh, MeshNamePar, ThisPe, NumPEs,IsParallel ) 
   
   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel

   ! Currently only one strategy to get the mesh is implemented 
   ! but there could be others.
   !
   ! This has not yet been tested in parallel and for sure
   ! it does not work for halo elements. 
   !-----------------------------------------------------------------
   CALL ElmerAsciiMesh( Step, PMesh, MeshNamePar, ThisPe, NumPEs, IsParallel ) 

 END SUBROUTINE LoadMeshStep

 !------------------------------------------------------------------------------
 ! Set the mesh dimension by studying the coordinate values.
 ! This could be less conservative also...
 !------------------------------------------------------------------------------    
 SUBROUTINE SetMeshDimension( Mesh )
   TYPE(Mesh_t), POINTER :: Mesh
   
   REAL(KIND=dp) :: x, y, z
   LOGICAL :: C(3)
   INTEGER :: i
   
   IF( Mesh % NumberOfNodes == 0 ) RETURN

   ! Compare value to some node, why not the 1st one
   x = Mesh % Nodes % x(1)
   y = Mesh % Nodes % y(1)
   z = Mesh % Nodes % z(1)
   
   C(1) = ANY( Mesh % Nodes % x /= x ) 
   C(2) = ANY( Mesh % Nodes % y /= y )  
   C(3) = ANY( Mesh % Nodes % z /= z )  

   ! This version is perhaps too liberal 
   Mesh % MeshDim = COUNT( C )
   Mesh % MaxDim = 0
   DO i=1,3
     IF( C(i) ) Mesh % MaxDim = i
   END DO
      
   CALL Info('SetMeshDimension','Dimension of mesh is: '//TRIM(I2S(Mesh % MeshDim)),Level=8)
   CALL Info('SetMeshDimension','Max dimension of mesh is: '//TRIM(I2S(Mesh % MaxDim)),Level=8)

 END SUBROUTINE SetMeshDimension

 
 !------------------------------------------------------------------------------
 !> Function to load mesh from disk.
 !------------------------------------------------------------------------------
 FUNCTION LoadMesh2( Model, MeshDirPar, MeshNamePar,&
     BoundariesOnly, NumProcs, MyPE, Def_Dofs, mySolver, &
     LoadOnly ) RESULT( Mesh )
   !------------------------------------------------------------------------------
   USE PElementMaps, ONLY : GetRefPElementNodes

   IMPLICIT NONE

   CHARACTER(LEN=*) :: MeshDirPar,MeshNamePar
   LOGICAL :: BoundariesOnly    
   INTEGER, OPTIONAL :: numprocs,mype,Def_Dofs(:,:), mySolver
   TYPE(Mesh_t),  POINTER :: Mesh
   TYPE(Model_t) :: Model
   LOGICAL, OPTIONAL :: LoadOnly 
   !------------------------------------------------------------------------------    
   INTEGER :: i,j,k,n
   INTEGER :: BaseNameLen, Save_Dim
   LOGICAL :: GotIt, Found
   CHARACTER(MAX_NAME_LEN) :: FileName
   TYPE(Element_t), POINTER :: Element
   TYPE(Matrix_t), POINTER :: Projector
   LOGICAL :: parallel, LoadNewMesh
   CHARACTER(*), PARAMETER :: Caller='LoadMesh'
   TYPE(ValueList_t), POINTER :: VList

   Mesh => Null()
   
   n = LEN_TRIM(MeshNamePar)
   DO WHILE (MeshNamePar(n:n)==CHAR(0).OR.MeshNamePar(n:n)==' ')
     n=n-1
   END DO
   IF(NumProcs<=1) THEN
     INQUIRE( FILE=MeshNamePar(1:n)//'/mesh.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Fatal(Caller,'Requested mesh > '//MeshNamePar(1:n)//' < does not exist!')
     END IF
     CALL Info(Caller,'Loading serial mesh!',Level=8)
    
   ELSE
     INQUIRE( FILE=MeshNamePar(1:n)//'/partitioning.'// & 
         TRIM(i2s(Numprocs))//'/part.1.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Warn(Caller,'Requested mesh > '//MeshNamePar(1:n)//' < in partition '&
           //TRIM(I2S(Numprocs))//' does not exist!')
       RETURN
     END IF
     CALL Info(Caller,'Loading parallel mesh for '//TRIM(I2S(Numprocs))//' partitions',Level=8)
   END IF
     
   Parallel = .FALSE.
   IF ( PRESENT(numprocs) .AND. PRESENT(mype) ) THEN
     IF ( numprocs > 1 ) Parallel = .TRUE.
   END IF

   Mesh => AllocateMesh()

   ! Get sizes of mesh structures for allocation
   !--------------------------------------------------------------------
   CALL LoadMeshStep( 1, Mesh, MeshNamePar, mype, numprocs, Parallel )

   ! Initialize and allocate mesh structures
   !---------------------------------------------------------------------
   IF( BoundariesOnly ) Mesh % NumberOfBulkElements = 0
   CALL InitializeMesh( Mesh )

   ! Get the (x,y,z) coordinates
   !--------------------------------------------------------------------------
   CALL LoadMeshStep( 2 )
   ! Permute and scale the coordinates.
   ! This also finds the mesh dimension. It is needed prior to getting the 
   ! elementtypes since wrong permutation or dimension may spoil that. 
   !-------------------------------------------------------------------
   CALL MapCoordinates()
   
   ! Get the bulk elements: element types, body index, topology
   !--------------------------------------------------------------------------
   CALL LoadMeshStep( 3 )

   ! Get the boundary elements: boundary types, boundary index, parents, topology
   !------------------------------------------------------------------------------
   CALL LoadMeshStep( 4 )

   ! Read elemental data - this is rarely used, parallel implementation lacking?
   !--------------------------------------------------------------------------
   i = LEN_TRIM(MeshNamePar)
   DO WHILE(MeshNamePar(i:i) == CHAR(0))
     i=i-1
   END DO
   BaseNameLen = i
   
   FileName = MeshNamePar(1:BaseNameLen)//'/mesh.elements.data'
   CALL ReadElementPropertyFile( FileName, Mesh )

   ! Read mesh.names - this could be saved by some mesh formats
   !--------------------------------------------------------------------------
   IF( ListGetLogical( Model % Simulation,'Use Mesh Names',Found ) ) THEN
     FileName = MeshNamePar(1:BaseNameLen)//'/mesh.names'
     CALL ReadTargetNames( Model, FileName )
   END IF


   ! Map bodies using Target Bodies and boundaries using Target Boundaries.
   ! This must be done before the element definitions are studied since
   ! then the pointer should be to the correct body index. 
   !------------------------------------------------------------------------
   CALL MapBodiesAndBCs()

   ! Read parallel mesh information: shared nodes
   !------------------------------------------------------------------
   CALL LoadMeshStep( 5 )

   ! Create the discontinuous mesh that accounts for the jumps in BCs
   ! This must be created after the whole mesh has been read in and 
   ! bodies and bcs have been mapped to full operation.
   ! To consider non-nodal elements it must be done before them.
   !--------------------------------------------------------------------
   CALL CreateDiscontMesh(Model,Mesh)

   !CALL CreateIntersectionBCs(Model,Mesh)
  
   ! Deallocate some stuff no longer needed
   !------------------------------------------------------------------
   CALL LoadMeshStep( 6 )

   CALL Info(Caller,'Loading mesh done',Level=8)
   
   IF( PRESENT( LoadOnly ) ) THEN
     CALL Info(Caller,'Only loading mesh, saving final preparation for later!',Level=12)     
     IF( LoadOnly ) RETURN
   END IF

   IF( PRESENT( mySolver ) ) THEN     
     VList => Model % Solvers(mySolver) % Values
   ELSE
     VList => Model % Simulation
   END IF
   IF(.NOT. ListGetLogical( VList,'Finalize Meshes Before Extrusion',Found ) ) THEN
     ! The final preparation for the mesh (including dof definitions) will be
     ! done only after the mesh has been extruded. 
     IF( ListCheckPresent( VList,'Extruded Mesh Levels') .OR. &
       ListCheckPresent( VList,'Extruded Mesh Layers') ) THEN
       CALL Info(Caller,'This mesh will be extruded, skipping finalization',Level=12)
       RETURN
     END IF
   END IF
     
   ! Prepare the mesh for next steps.
   ! For example, create non-nodal mesh structures, periodic projectors etc. 
   CALL PrepareMesh(Model,Mesh,Parallel,Def_Dofs,mySolver)      
   CALL Info(Caller,'Preparing mesh done',Level=8)

   
 CONTAINS


   !------------------------------------------------------------------------------
   ! Map bodies and boundaries as prescirbed by the 'Target Bodies' and 
   ! 'Target Boundaries' keywords.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapBodiesAndBCs()

     TYPE(Element_t), POINTER :: Element
     INTEGER, ALLOCATABLE :: IndexMap(:), TmpIndexMap(:)
     INTEGER, POINTER :: Blist(:)
     INTEGER :: id,minid,maxid,body,bndry,DefaultTargetBC, DefaultTargetBody


     ! If "target bodies" is used map the bodies accordingly
     !------------------------------------------------------
     Found = .FALSE. 
     DefaultTargetBody = 0
     DO id=1,Model % NumberOfBodies
       IF( ListCheckPresent( Model % Bodies(id) % Values,'Target Bodies') ) Found = .TRUE.
       IF(ListGetLogical( Model % Bodies(id) % Values, &
           'Default Target', GotIt)) THEN
         DefaultTargetBody = id
         Found = .TRUE.
       END IF
     END DO

     IF( DefaultTargetBody /= 0 ) THEN
       CALL Info('MapBodiesAndBCs','Default Target Body: '&
           //TRIM(I2S(DefaultTargetBody)),Level=8)
     END IF
     
     IF( Found ) THEN
       CALL Info('MapBodiesAndBCs','Remapping bodies',Level=8)      
       minid = HUGE( minid ) 
       maxid = -HUGE( maxid ) 
       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         id = Element % BodyId
         minid = MIN( id, minid ) 
         maxid = MAX( id, maxid )
       END DO
       IF( minid > maxid ) THEN
         CALL Fatal('MapBodiesAndBCs','Body indexes are screwed!')
       END IF
       CALL Info('MapBodiesAndBCs','Minimum initial body index: '//TRIM(I2S(minid)),Level=6 )
       CALL Info('MapBodiesAndBCs','Maximum initial body index: '//TRIM(I2S(maxid)),Level=6 )

       minid = MIN( 1, minid ) 
       maxid = MAX( Model % NumberOfBodies, maxid ) 
       ALLOCATE( IndexMap(minid:maxid) )
       IndexMap = 0

       DO id=1,Model % NumberOfBodies
         BList => ListGetIntegerArray( Model % Bodies(id) % Values, &
             'Target Bodies', GotIt ) 
         IF ( Gotit ) THEN
           DO k=1,SIZE(BList)
             body = Blist(k)
             IF( body > maxid .OR. body < minid ) THEN
#if 0
               CALL Warn('MapBodiesAndBCs','Unused body entry in > Target Bodies <  : '&
                   //TRIM(I2S(body)) )              
#endif
             ELSE IF( IndexMap( body ) /= 0 ) THEN
               CALL Warn('MapBodiesAndBCs','Multiple bodies have same > Target Bodies < entry : '&
                   //TRIM(I2S(body)))
             ELSE
               IndexMap( body ) = id 
             END IF
           END DO
         ELSE
           IF( DefaultTargetBody == 0 ) THEN
             IF( IndexMap( id ) /= 0 ) THEN
               CALL Warn('MapBodiesAndBCs','Unset body already set by > Target Boundaries < : '&
                   //TRIM(I2S(id)) )
             ELSE 
               IndexMap( id ) = id
             END IF
           END IF
         END IF
           
       END DO

       IF( .FALSE. ) THEN
         PRINT *,'Body mapping'
         DO id=minid,maxid
           IF( IndexMap( id ) /= 0 ) PRINT *,id,' : ',IndexMap(id)
         END DO
       END IF

       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         id = Element % BodyId

         IF( IndexMap( id ) == 0 ) THEN
           IF( DefaultTargetBody /= 0 ) THEN
             IndexMap( id ) = DefaultTargetBody
           END IF
         END IF

         Element % BodyId = IndexMap( id ) 
       END DO

       DEALLOCATE( IndexMap )
     ELSE
       CALL Info('MapBodiesAndBCs','Skipping remapping of bodies',Level=10)      
     END IF


     IF( Mesh % NumberOfBoundaryElements == 0 ) RETURN

     ! Target boundaries are usually given so this is not conditional
     !---------------------------------------------------------------
     CALL Info('MapBodiesAndBCs','Remapping boundaries',Level=8)      
     minid = HUGE( minid ) 
     maxid = -HUGE( maxid ) 
     DO i=Mesh % NumberOfBulkElements+1,&
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       id = Element % BoundaryInfo % Constraint
       minid = MIN( id, minid ) 
       maxid = MAX( id, maxid )
     END DO


     CALL Info('MapBodiesAndBCs','Minimum initial boundary index: '//TRIM(I2S(minid)),Level=6 )
     CALL Info('MapBodiesAndBCs','Maximum initial boundary index: '//TRIM(I2S(maxid)),Level=6 )
     IF( minid > maxid ) THEN
       CALL Fatal('MapBodiesAndBCs','Boundary indexes are screwed')
     END IF

     minid = MIN( minid, 1 ) 
     maxid = MAX( maxid, Model % NumberOfBCs ) 
     ALLOCATE( IndexMap(minid:maxid) )
     IndexMap = 0


     DO j=1,Model % NumberOfBoundaries
       id = ListGetInteger( Model % Boundaries(j) % Values, &
           'Boundary Condition',GotIt, minv=1, maxv=Model % NumberOFBCs )
       IF( id == 0 ) CYCLE
       bndry = Model % BoundaryId(j)
       IF( bndry > maxid ) THEN
         CALL Warn('MapBodiesAndBCs','BoundaryId exceeds range')
       ELSE IF( bndry == 0 ) THEN
         CALL Warn('MapBodiesAndBCs','BoundaryId is zero')
       ELSE
         IndexMap( bndry ) = id
       END IF
     END DO

     DefaultTargetBC = 0
     DO id=1,Model % NumberOfBCs
       IF(ListGetLogical( Model % BCs(id) % Values, &
           'Default Target', GotIt)) DefaultTargetBC = id       
       BList => ListGetIntegerArray( Model % BCs(id) % Values, &
           'Target Boundaries', GotIt )
       IF ( Gotit ) THEN
         DO k=1,SIZE(BList)
           bndry = Blist(k)
           IF( bndry > maxid ) THEN
#if 0
  in my opinion, this is quite usual ... Juha
             CALL Warn('MapBodiesAndBCs','Unused BC entry in > Target Boundaries <  : '&
                 //TRIM(I2S(bndry)) )              
#endif
           ELSE IF( IndexMap( bndry ) /= 0 ) THEN
             CALL Warn('MapBodiesAndBCs','Multiple BCs have same > Target Boundaries < entry : '&
                 //TRIM(I2S(bndry)) )
           ELSE 
             IndexMap( bndry ) = id 
           END IF
         END DO
       ELSE
         IF (ListCheckPresent(Model % BCs(id) % Values, 'Target Nodes') .OR. &
             ListCheckPresent(Model % BCs(id) % Values, 'Target Coordinates')) &
             CYCLE
         IF (IndexMap( id ) /= 0 .AND. id == DefaultTargetBC ) THEN ! DefaultTarget has been given
           CALL Warn('MapBodiesAndBCs','Default Target is a Target Boundaries entry in > Boundary Condition < : '&
               //TRIM(I2S(IndexMap(id))) )
         END IF
         !
         !IF( IndexMap( id ) /= 0 .AND. id /= DefaultTargetBC ) THEN
         !  CALL Warn(Caller,'Unset BC already set by > Target Boundaries < : '&
         !      //TRIM(I2S(id)) )
         !ELSE 
         !  ! IndexMap( id ) = id
         !END IF
       END IF
     END DO

     IF( .FALSE. ) THEN
       PRINT *,'Boundary mapping'
       DO id=minid,maxid
         IF( IndexMap( id ) /= 0 ) PRINT *,id,' : ',IndexMap(id)
       END DO
     END IF

     IF( DefaultTargetBC /= 0 ) THEN
       CALL Info('MapBodiesAndBCs','Default Target BC: '&
           //TRIM(I2S(DefaultTargetBC)),Level=8)
     END IF


     DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       n = Element % TYPE % NumberOfNodes
       bndry = Element % BoundaryInfo % Constraint 

       IF( bndry > maxid .OR. bndry < minid ) THEN
         CALL Warn('MapBodiesAndBCs','Boundary index '//TRIM(I2S(bndry))&
             //' not in range: '//TRIM(I2S(minid))//','//TRIM(I2S(maxid)) )
       END IF

       IF( IndexMap( bndry ) < 0 ) THEN
         Element % BoundaryInfo % Constraint = 0
         CYCLE

       ELSE IF( IndexMap( bndry ) == 0 ) THEN
         IF( DefaultTargetBC /= 0 ) THEN
!          PRINT *,'Default boundary map: ',bndry,DefaultTargetBC
           IndexMap( bndry ) = DefaultTargetBC
         ELSE 
!          IF( bndry <= Model % NumberOfBCs ) THEN            
!            PRINT *,'Unmapped boundary: ',bndry
!          ELSE
!            PRINT *,'Unused boundary: ',bndry
!          END IF
           IndexMap( bndry ) = -1 
           Element % BoundaryInfo % Constraint = 0           
           CYCLE
         END IF
       END IF

       bndry = IndexMap( bndry ) 
       Element % BoundaryInfo % Constraint = bndry 

       IF( bndry <= Model % NumberOfBCs ) THEN
         Element % BodyId  = ListGetInteger( &
             Model % BCs(bndry) % Values, 'Body Id', Gotit, 1, Model % NumberOfBodies )
         Element % BoundaryInfo % OutBody = &
             ListGetInteger( Model % BCs(bndry) % Values, &
             'Normal Target Body', GotIt, maxv=Model % NumberOFBodies ) 
       END IF
     END DO

     DEALLOCATE( IndexMap ) 

   END SUBROUTINE MapBodiesAndBCs

   

   !------------------------------------------------------------------------------
   ! Map and scale coordinates, and increase the size of the coordinate
   ! vectors, if requested.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapCoordinates()

     REAL(KIND=dp), POINTER CONTIG :: NodesX(:), NodesY(:), NodesZ(:)
     REAL(KIND=dp), POINTER :: Wrk(:,:)
     INTEGER, POINTER :: CoordMap(:)
     REAL(KIND=dp) :: CoordScale(3)
     INTEGER :: mesh_dim, model_dim
     
     ! Perform coordinate mapping
     !------------------------------------------------------------
     CoordMap => ListGetIntegerArray( Model % Simulation, &
         'Coordinate Mapping',GotIt )
     IF ( GotIt ) THEN
       CALL Info('MapCoordinates','Performing coordinate mapping',Level=8)

       IF ( SIZE( CoordMap ) /= 3 ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'MapCoordinates', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'MapCoordinates', Message )
       END IF

       IF ( ALL( CoordMap(1:3) /= 1 ) .OR. ALL( CoordMap(1:3) /= 2 ) .OR. ALL( CoordMap(1:3) /= 3 ) ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'MapCoordinates', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'MapCoordinates', Message )
       END IF

       IF( CoordMap(1) == 1 ) THEN
         NodesX => Mesh % Nodes % x
       ELSE IF( CoordMap(1) == 2 ) THEN
         NodesX => Mesh % Nodes % y
       ELSE
         NodesX => Mesh % Nodes % z
       END IF

       IF( CoordMap(2) == 1 ) THEN
         NodesY => Mesh % Nodes % x
       ELSE IF( CoordMap(2) == 2 ) THEN
         NodesY => Mesh % Nodes % y
       ELSE
         NodesY => Mesh % Nodes % z
       END IF

       IF( CoordMap(3) == 1 ) THEN
         NodesZ => Mesh % Nodes % x
       ELSE IF( CoordMap(3) == 2 ) THEN
         NodesZ => Mesh % Nodes % y
       ELSE
         NodesZ => Mesh % Nodes % z
       END IF

       Mesh % Nodes % x => NodesX
       Mesh % Nodes % y => NodesY
       Mesh % Nodes % z => NodesZ
     END IF

     ! Determine the mesh dimension 
     !----------------------------------------------------------------------------
     CALL SetMeshDimension( Mesh )
     
     mesh_dim = Mesh % MaxDim

     ! Scaling of coordinates
     !-----------------------------------------------------------------------------
     Wrk => ListGetConstRealArray( Model % Simulation,'Coordinate Scaling',GotIt )    
     IF( GotIt ) THEN            
       CoordScale = 1.0_dp
       DO i=1,mesh_dim
         j = MIN( i, SIZE(Wrk,1) )
         CoordScale(i) = Wrk(j,1)
       END DO
       WRITE(Message,'(A,3ES10.3)') 'Scaling coordinates:',CoordScale(1:3)
       CALL Info('MapCoordinates',Message) 
       Mesh % Nodes % x = CoordScale(1) * Mesh % Nodes % x
       IF( mesh_dim > 1 ) Mesh % Nodes % y = CoordScale(2) * Mesh % Nodes % y
       IF( mesh_dim > 2 ) Mesh % Nodes % z = CoordScale(3) * Mesh % Nodes % z
     END IF

   END SUBROUTINE MapCoordinates

 !------------------------------------------------------------------------------
 END FUNCTION LoadMesh2
 !------------------------------------------------------------------------------


 !> Prepare a clean nodal mesh as it comes after being loaded from disk.
 !> Study the non-nodal elements (face, edge, DG, and p-elements)
 !> Create parallel info for the non-nodal elements
 !> Enlarge the coordinate vectors for p-elements.
 !> Generate static projector for periodic BCS.
 !-------------------------------------------------------------------
 SUBROUTINE PrepareMesh( Model, Mesh, Parallel, Def_Dofs, mySolver )

   TYPE(Model_t) :: Model
   TYPE(Mesh_t), POINTER :: Mesh
   LOGICAL :: Parallel
   INTEGER, OPTIONAL :: Def_Dofs(:,:), mySolver
   LOGICAL :: Found
   CHARACTER(*),PARAMETER :: Caller='PrepareMesh'

      
   IF( Mesh % MaxDim == 0) THEN
     CALL SetMeshDimension( Mesh )
   END IF
   Model % DIMENSION = MAX( Model % DIMENSION, Mesh % MaxDim ) 
     
   CALL CreateIntersectionBCs(Model,Mesh)

   CALL NonNodalElements()

   IF( Parallel ) THEN
     CALL ParallelNonNodalElements()
   END IF
     
   CALL EnlargeCoordinates( Mesh ) 

   CALL FollowCurvedBoundary( Model, Mesh, .FALSE. ) 
     
   CALL GeneratePeriodicProjectors( Model, Mesh )    
   
   IF( ListGetLogical( Model % Simulation,'Inspect Quadratic Mesh', Found ) ) THEN
     CALL InspectQuadraticMesh( Mesh ) 
   END IF
   
   IF( ListGetLogical( Model % Simulation,'Inspect Mesh',Found ) ) THEN
     CALL InspectMesh( Mesh ) 
   END IF

   IF(ListGetLogical( Model % Simulation, 'Parallel Reduce Element Max Sizes', Found ) ) THEN
     Mesh % MaxElementDOFs  = ParallelReduction( Mesh % MaxElementDOFs,2  ) 
     Mesh % MaxElementNodes = ParallelReduction( Mesh % MaxElementNodes,2 ) 
   END IF

   
   
 CONTAINS
     

   ! Check for the non-nodal element basis
   !--------------------------------------------------------
   SUBROUTINE NonNodalElements()

     INTEGER, POINTER :: EdgeDofs(:), FaceDofs(:)
     INTEGER :: i, j, k, l, s, n, DGIndex, body_id, body_id0, eq_id, solver_id, el_id, &
         mat_id
     LOGICAL :: NeedEdges, Found, FoundDef0, FoundDef, FoundEq, GotIt, MeshDeps, &
         FoundEqDefs, FoundSolverDefs(Model % NumberOfSolvers), &
         FirstOrderElements, InheritDG, Hit, Stat, &
         UpdateDefDofs(Model % NumberOfSolvers)
     TYPE(Element_t), POINTER :: Element, Parent, pParent
     TYPE(Element_t) :: DummyElement
     TYPE(ValueList_t), POINTER :: Vlist
     INTEGER :: inDOFs(10,6)
     CHARACTER(MAX_NAME_LEN) :: ElementDef0, ElementDef
     
     
     EdgeDOFs => NULL()
     CALL AllocateVector( EdgeDOFs, Mesh % NumberOfBulkElements, Caller )
     FaceDOFs => NULL()
     CALL AllocateVector( FaceDOFs, Mesh % NumberOfBulkElements, Caller )     
    
     DGIndex = 0

     InDofs = 0
     InDofs(:,1) = 1
     IF ( PRESENT(Def_Dofs) ) THEN
       inDofs = Def_Dofs
     ELSE
       DO s=1,Model % NumberOfSolvers
         DO i=1,6
           DO j=1,10
             inDofs(j,i) = MAX(Indofs(j,i),MAXVAL(Model % Solvers(s) % Def_Dofs(j,:,i)))
           END DO
         END DO
       END DO
     END IF

     ! P-basis only over 1st order elements:
     ! -------------------------------------
     FirstOrderElements = .TRUE.
     DO i=1,Mesh % NumberOfBulkElements
       IF (Mesh % Elements(i) % Type % BasisFunctionDegree>1) THEN
         FirstOrderElements = .FALSE.; EXIT
       END IF
     END DO

    !
    ! Check whether the "Element" definitions can depend on mesh
    ! -----------------------------------------------------------
    MeshDeps = .FALSE.  ! The order of p-basis given with a MATC function
    FoundEqDefs = .FALSE.;  FoundSolverDefs = .FALSE.

    !
    ! As a preliminary step, check if an element definition is given 
    ! in an equation section. The more common way is to give the element
    ! definition in a solver section.
    !
    DO eq_id=1,Model % NumberOFEquations
      Vlist => Model % Equations(eq_id) % Values
      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0)
      FoundEqDefs = FoundEqDefs .OR. FoundDef0

      IF (FoundDef0) THEN
        !
        ! Check if the order of p-basis is defined by calling a special
        ! MATC function:
        !
        j = INDEX(ElementDef0,'p:')
        IF (j>0 .AND. ElementDef0(j+2:j+2)=='%') MeshDeps = .TRUE.
      ELSE
        !
        ! Check if element definitions are given for each solver separately
        ! by using a special keyword construct and tag the corresponding
        ! entries in the list of the solvers. 
        ! 
        DO Solver_id=1,Model % NumberOfSolvers
          IF (PRESENT(mySolver)) THEN
            IF ( Solver_id /= mySolver ) CYCLE
          ELSE
            ! Respect definitions given in the solver section:
            IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
          END IF

          ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)
          FoundSolverDefs(Solver_id) = FoundSolverDefs(solver_id) .OR. FoundDef

          IF (FoundDef) THEN
            j = INDEX(ElementDef,'p:')
            IF (j>0 .AND. ElementDef(j+2:j+2)=='%') MeshDeps = .TRUE.
          END IF
        END DO
      END IF
    END DO

    !
    ! Tag solvers for which the element definition has been given in
    ! a solver section. The function LoadModel has already read these
    ! element definitions except for cases where the order of p-basis is
    ! defined in terms of a MATC function. The array UpdateDefDofs will
    ! show whether element definitions should be re-read.
    !
    UpdateDefDofs = .TRUE.
    DO solver_id=1,Model % NumberOfSolvers
      Vlist => Model % Solvers(solver_id) % Values

      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0)

      IF (FoundDef0) THEN
        FoundSolverDefs(Solver_id) = .TRUE.

        j = INDEX(ElementDef0,'p:')
        IF (j>0 .AND. ElementDef0(j+2:j+2)=='%') THEN
          meshdeps = .TRUE.
        ELSE
          ! Solverwise element definitions have already be read in LoadModel,
          ! indicate that re-reading is not needed here
          UpdateDefDofs(Solver_id) = .FALSE.
        END IF
      END IF
    END DO

    ! The basic case without the order of p-basis being defined by a MATC function:
    !
    IF (.NOT.MeshDeps) THEN
      FoundDef0 = .FALSE.
      DO body_id=1,Model % NumberOfBodies
        ElementDef0 = ' '
        Vlist => Model % Bodies(body_id) % Values
        eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
        IF ( FoundEq ) THEN
          Vlist => Model % Equations(eq_id) % Values
          IF (FoundEqDefs) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

          DO solver_id=1,Model % NumberOfSolvers

            IF(PRESENT(mySolver)) THEN
              IF ( Solver_id /= mySolver ) CYCLE
            ELSE
              IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
            END IF

            FoundDef = .FALSE.
            IF(FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)

            IF ( FoundDef ) THEN
              CALL GetMaxDefs( Model, Mesh, DummyElement, ElementDef, solver_id, body_id, Indofs )
            ELSE
              IF (UpdateDefDofs(Solver_id)) THEN
                IF (.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) &
                    ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

                CALL GetMaxDefs( Model, Mesh, DummyElement, ElementDef0, solver_id, body_id, Indofs )

                IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
              ! ELSE
              !   PRINT *, 'NO NEED TO RECREATE DEF_DOFS '
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF

     ! non-nodal elements in bulk elements
     !------------------------------------------------------------
     body_id0 = -1; FoundDef=.FALSE.; FoundEq=.FALSE.
     ElementDef = ' '

     !
     ! Check whether face DOFs have been generated by "-quad_face b: ..." or
     ! "-tri_face b: ..."
     !
     NeedEdges = ANY( inDOFs(9:10,5)>0 )

     DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       body_id = Element % BodyId
       n = Element % TYPE % NumberOfNodes
       
       ! Check if the order of p-basis depends on a MATC function
       IF ( Meshdeps ) THEN
         IF ( body_id/=body_id0 ) THEN
           Vlist => Model % Bodies(body_id) % Values
           eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
           ElementDef0 = ' '
         END IF

         IF ( FoundEq ) THEN
           Vlist => Model % Equations(eq_id) % Values
           FoundDef0 = .FALSE.
           IF( FoundEqDefs.AND.body_id/=body_id0 ) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

           DO solver_id=1,Model % NumberOfSolvers
             IF(PRESENT(mySolver)) THEN
               IF ( Solver_id /= mySolver ) CYCLE
             ELSE
               IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
             END IF

             FoundDef = .FALSE.
             IF (FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)

             IF ( FoundDef ) THEN
               CALL GetMaxDefs( Model, Mesh, Element, ElementDef, solver_id, body_id, Indofs )
             ELSE
               IF (UpdateDefDofs(Solver_id)) THEN
                 IF (.NOT. FoundDef0.AND.FoundSolverDefs(solver_id)) &
                     ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

                 CALL GetMaxDefs( Model, Mesh, Element, ElementDef0, solver_id, body_id, Indofs )

                 IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
               END IF
             END IF
           END DO
         END IF
         body_id0 = body_id
       END IF


       el_id = Element % TYPE % ElementCode / 100

       ! Apply the elementtypes

       Element % NDOFs = n * MAX(0,inDOFs(el_id,1)) ! The count of all nodal DOFs for the element
       EdgeDOFs(i) = MAX(0,inDOFs(el_id,2))
       FaceDOFs(i) = MAX(0,inDOFs(el_id,3))

       IF ( inDofs(el_id,4) == 0 ) THEN
         inDOFs(el_id,4) = n
       END IF

       NULLIFY( Element % DGIndexes )
       IF ( inDOFs(el_id,4) > 0 ) THEN
         CALL AllocateVector( Element % DGIndexes, inDOFs(el_id,4))
         IF( indofs(el_id,4) /= Element % TYPE % NumberOfNodes ) &
             PRINT *,'Element:',Element % TYPE % ElementCode, indofs(el_id,4)
         DO j=1,inDOFs(el_id,4)
           DGIndex = DGIndex + 1
           Element % DGIndexes(j) = DGIndex
         END DO
       END IF
       Element % DGDOFs = MAX(0,inDOFs(el_id,4))
       NeedEdges = NeedEdges .OR. ANY( inDOFs(el_id,2:4)>0 )
       
       
       ! Check if given element is a p element
       IF (FirstOrderElements .AND. inDOFs(el_id,6) > 0) THEN
         CALL AllocatePDefinitions(Element)
         NeedEdges = .TRUE.
         
         ! Calculate element bubble dofs and set element p

         Element % PDefs % P = inDOFs(el_id,6)   ! NOTE: If the order of p-basis is given by
                                                 ! a MATC function, the order is here defined
                                                 ! to be the maximum order over the element
                                                 ! processed so far. This is 
                                                 ! erroneous as the resulting p-distribution  
                                                 ! thus depends on the numbering of geometric
                                                 ! entities.
         !
         ! Try to fix the issue described in the above remark in a special case 
         ! where a single element definition is given in the equation section:
         !
         IF (FoundEqDefs .AND. Model % NumberOfSolvers > 0) THEN
           ! All solvers have the same element definition, pick one of these
           ! to set the polynomial degree:
           Element % PDefs % P = Model % Solvers(1) % Def_Dofs(el_id,Body_Id,6)
         END IF

         IF ( inDOFs(el_id,5) > 0 ) THEN
           Element % BDOFs = inDOFs(el_id,5)
         ELSE
           Element % BDOFs = getBubbleDOFs(Element, Element % PDefs % P)
         END IF

         ! All elements in actual mesh are not edges
         Element % PDefs % pyramidQuadEdge = .FALSE.
         Element % PDefs % isEdge = .FALSE.

         ! If element is of type tetrahedron and is a p element, 
         ! do the Ainsworth & Coyle trick
         IF (Element % TYPE % ElementCode == 504) CALL ConvertToACTetra(Element)
         CALL GetRefPElementNodes( Element % Type,  Element % Type % NodeU, &
             Element % Type % NodeV, Element % Type % NodeW )
       ELSE 
         ! Clear P element definitions and set manual bubbles
         Element % PDefs => NULL()
         Element % BDOFs = MAX(0,inDOFs(el_id,5))
         ! WRITE (*,*) Element % BDOFs
       END IF

       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes,Element % TYPE % NumberOfNodes )
     END DO

     InheritDG = .FALSE.
     IF( dgindex > 0 ) THEN
       InheritDG = ListCheckPresentAnyMaterial( CurrentModel,'DG Parent Material')
     END IF
     
     ! non-nodal elements in boundary elements
     !------------------------------------------------------------    
     DO i = Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('NonNodalElements','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
         CALL Fatal('NonNodalElements','Type in Element '//TRIM(I2S(i))//' not associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       el_id = Element % TYPE % ElementCode / 100
       Element % NDOFs  = n * MAX(0,inDOFs(el_id,1))
       
       !
       ! NOTE: The following depends on what dofs have been introduced
       ! by using the construct "-quad_face b: ..." and
       ! "-tri_face b: ..."
       !
       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) THEN
         IF( Element % BoundaryInfo % Left % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           Element % BDOFs = &
               EdgeDOFs(Element % BoundaryInfo % Left % ElementIndex)
         ELSE
           Element % BDOFs = FaceDOFs(Element % BoundaryInfo % Left % ElementIndex)
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) THEN
         IF ( Element % BoundaryInfo % Right % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           Element % BDOFs = &
               EdgeDOFs(Element % BoundaryInfo % Right % ElementIndex)
         ELSE
           Element % BDOFs = FaceDOFs(Element % BoundaryInfo % Right % ElementIndex)
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF

       ! Optionally also set DG indexes for BCs
       ! It is easy for outside boundaries, but for internal boundaries
       ! we need a flag "DG Parent Material".
       IF( InheritDG ) THEN
         IF(.NOT. ASSOCIATED( Element % DGIndexes ) ) THEN
           ALLOCATE( Element % DGIndexes(n) )
           Element % DGIndexes = 0
         END IF
         
         Hit = .TRUE.
         k = 0
         DO l=1,2        
           IF(l==1) THEN
             Parent => Element % BoundaryInfo % Left
           ELSE
             Parent => Element % BoundaryInfo % Right
           END IF
           IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
           k = k + 1
           pParent => Parent
           
           mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,&
               'Material',Found )
           IF(mat_id > 0 ) THEN           
             VList => CurrentModel % Materials(mat_id) % Values
           END IF
           IF( ASSOCIATED(Vlist) ) THEN
             Hit = ListGetLogical(Vlist,'DG Parent Material',Found )
           END IF
           IF( Hit ) EXIT
         END DO
         
         IF( k == 0 ) THEN
           CALL Fatal('NonnodalElements','Cannot define DG indexes for BC!')
         ELSE IF( k == 1 ) THEN
           Parent => pParent        
         ELSE IF(.NOT. Hit ) THEN
           CALL Fatal('NonnodalElements','Cannot define DG indexes for internal BC!')       
         END IF
         
         DO l=1,n
           DO j=1, Parent % TYPE % NumberOfNodes
             IF( Element % NodeIndexes(l) == Parent % NodeIndexes(j) ) THEN
               Element % DGIndexes(l) = Parent % DGIndexes(j)
               EXIT
             END IF
           END DO
         END DO
       END IF
       
     END DO

     IF ( Mesh % MaxElementDOFs <= 0 ) Mesh % MaxElementDOFs = Mesh % MaxElementNodes 

     ! Override automated "NeedEdges" if requested by the user.
     !------------------------------------------------------------------------------------
     IF(PRESENT(mySolver)) THEN
       Stat = ListGetLogical(Model % Solvers(mySolver) % Values, 'Need Edges', Found)
       IF(Found) NeedEdges = Stat

       IF( ListGetLogical(Model % Solvers(mySolver) % Values, 'NeedEdges', Found) ) THEN
         IF(.NOT. NeedEdges) CALL Fatal('NonNodalElements','Use "Need Edges" instead of "NeedEdges"') 
       END IF
     END IF

     IF( Mesh % MeshDim == 2 ) THEN
       Stat = ListGetLogical(Model % Simulation, 'Need Edges 2D', Found)
       IF(Found) NeedEdges = Stat
     END IF

     IF( Mesh % MeshDim == 3 ) THEN
       Stat = ListGetLogical(Model % Simulation, 'Need Edges 3D', Found)
       IF(Found) NeedEdges = Stat
     END IF
     
     IF ( NeedEdges ) THEN
       CALL Info('NonNodalElements','Requested elements require creation of edges',Level=8)
       CALL SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs)
     END IF

     CALL SetMeshMaxDOFs(Mesh)

     IF( ASSOCIATED(EdgeDOFs) ) DEALLOCATE(EdgeDOFs )
     IF( ASSOCIATED(FaceDOFs) ) DEALLOCATE(FaceDOFs)

     IF( Mesh % MaxFaceDofs > 0 ) THEN
       CALL Info('NonNodalElements','Face dofs max: '//TRIM(I2S(Mesh % MaxFaceDofs)),Level=12)
     END IF
     IF( Mesh % MaxEdgeDofs > 0 ) THEN
       CALL Info('NonNodalElements','Edge dofs max: '//TRIM(I2S(Mesh % MaxEdgeDofs)),Level=12)
     END IF
     IF( Mesh % MaxElementDofs > 0 ) THEN
       CALL Info('NonNodalElements','Element dofs max: '//TRIM(I2S(Mesh % MaxElementDofs)),Level=12)
     END IF

   END SUBROUTINE NonNodalElements


   ! When the parallel nodal neighbours have been found 
   ! perform numbering for face and edge elements as well.
   !-------------------------------------------------------------------    
   SUBROUTINE ParallelNonNodalElements()

     INTEGER :: i,j,k,n     
     TYPE(Element_t), POINTER :: Element

     ! To be on the safe side create the parallel info if it is missing.
     IF( Mesh % NumberOfNodes > 0 ) THEN
       n = SIZE( Mesh % ParallelInfo % NeighbourList )              
       ! For unset neighbours just set the this partition to be the only owner
       DO i=1,n
         IF (.NOT.ASSOCIATED(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) THEN
           CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(i) % Neighbours,1)
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % mype
         END IF
       END DO
     END IF
       
     ! Create parallel numbering of faces
     CALL SParFaceNumbering(Mesh, .TRUE. )

     ! Create parallel numbering for edges
     CALL SParEdgeNumbering(Mesh, .TRUE.)

     ! There are mainly implemented for parallel debugging.
     ! The whole sequence is only activated when "Max Output Level >= 10". 
     IF( InfoActive(10) ) THEN
       CALL Info('ParallelNonNodalElements','Number of initial nodes: '&
           //TRIM(I2S(Mesh % NumberOfNodes)))
       
       CALL Info('ParallelNonNodalElements','Number of initial faces: '&
           //TRIM(I2S(Mesh % NumberOfFaces)))
       
       CALL Info('ParallelNonNodalElements','Number of initial edges: '&
           //TRIM(I2S(Mesh % NumberOfEdges)))
       
       j = 0; k = 0
       DO i=1,Mesh % NumberOfNodes
         IF( SIZE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) > 1 ) THEN
           j = j + 1
           IF( Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPe ) k = k + 1
         END IF
       END DO      
       CALL Info('ParallelNonNodalElements','Number of shared nodes: '//TRIM(I2S(j)))
       CALL Info('ParallelNonNodalElements','Number of owned shared nodes: '//TRIM(I2S(k)))
            
       IF( Mesh % NumberOfFaces > 0 ) THEN
         j = 0; k = 0 
         DO i=1,Mesh % NumberOfFaces
           IF( SIZE( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours ) > 1 ) THEN
             j = j + 1 
             IF( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours(1) == ParEnv % MyPe ) k = k + 1   
           END IF
         END DO
         CALL Info('ParallelNonNodalElements','Number of shared faces: '//TRIM(I2S(j)))
         CALL Info('ParallelNonNodalElements','Number of owned shared faces: '//TRIM(I2S(k)))

#if 0
         DO i=1,Mesh % NumberOfFaces
           IF( SIZE( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours ) == 1 ) THEN
             BLOCK
               TYPE(Element_t), POINTER :: Face
               Face => Mesh % Faces(i)
               k = 0
               DO j=1,Face % TYPE % NumberOfNodes 
                 IF( SIZE( Mesh % ParallelInfo % NeighbourList(Face % NodeIndexes(j)) % Neighbours ) > 1 ) k = k + 1 
               END DO
               IF( k == Face % TYPE % NumberOfNodes ) THEN
                 PRINT *,'Face is shared but not listed!',ParEnv % MyPe, Mesh % NumberOfFaces,i
               END IF
             END BLOCK
           ELSE
             PRINT *,'Face is shared and listed: ',ParEnv % MyPe, Mesh % NumberOfFaces,i             
           END IF
         END DO
#endif   

       END IF
       
       IF( Mesh % NumberOfEdges > 0 ) THEN
         j = 0; k = 0
         DO i=1,Mesh % NumberOfEdges
           IF( SIZE( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours ) > 1 ) THEN
             j = j + 1
             IF( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours(1) == ParEnv % MyPe ) k = k + 1   
           END IF
         END DO
         CALL Info('ParallelNonNodalElements','Number of shared edges: '//TRIM(I2S(j)))
         CALL Info('ParallelNonNodalElements','Number of owned shared edges: '//TRIM(I2S(k)))
       END IF
     END IF


     DO i=1,Mesh % NumberOfFaces
       Mesh % MinFaceDOFs = MIN(Mesh % MinFaceDOFs,Mesh % Faces(i) % BDOFs)
       Mesh % MaxFaceDOFs = MAX(Mesh % MaxFaceDOFs,Mesh % Faces(i) % BDOFs)
     END DO
     IF(Mesh % MinFaceDOFs > Mesh % MaxFaceDOFs) Mesh % MinFaceDOFs = Mesh % MaxFaceDOFs

     DO i=1,Mesh % NumberOfEdges
       Mesh % MinEdgeDOFs = MIN(Mesh % MinEdgeDOFs,Mesh % Edges(i) % BDOFs)
       Mesh % MaxEdgeDOFs = MAX(Mesh % MaxEdgeDOFs,Mesh % Edges(i) % BDOFs)
     END DO
     IF(Mesh % MinEdgeDOFs > Mesh % MaxEdgeDOFs) Mesh % MinEdgeDOFs = Mesh % MaxEdgeDOFs

     ! Set max element dofs here (because element size may have changed
     ! when edges and faces have been set). This is the absolute worst case.
     ! Element which has MaxElementDOFs may not even be present as a 
     ! real element
     DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)        
       Mesh % MaxElementDOFs = MAX( Mesh % MaxElementDOFs, &
           Element % TYPE % NumberOfNodes + &
           Element % TYPE % NumberOfEdges * Mesh % MaxEdgeDOFs + &
           Element % TYPE % NumberOfFaces * Mesh % MaxFaceDOFs + &
           Element % BDOFs, &
           Element % DGDOFs )
     END DO

   END SUBROUTINE ParallelNonNodalElements

   
 END SUBROUTINE PrepareMesh



  !-------------------------------------------------------------------------------
  !> Communicate logical tag related to mesh or linear system.
  !> This could related to setting Neumann BCs to zero, for example.
  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateParallelSystemTag(ParallelInfo,Ltag,Itag,CommVal)
  !-------------------------------------------------------------------------------
     TYPE (ParallelInfo_t), POINTER :: ParallelInfo
     LOGICAL, POINTER, OPTIONAL :: LTag(:)   !< Logical tag, if used
     INTEGER, POINTER, OPTIONAL :: ITag(:)   !< Integer tag, if used
     LOGICAL, OPTIONAL :: CommVal            !< If integer tag is used, should we consider also the value

     LOGICAL, POINTER :: IsNeighbour(:)
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:), s_i(:,:), r_i(:)
     INTEGER :: i,j,k,l,n,nn,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)
     INTEGER :: NewZeros, nsize
     LOGICAL :: UseL, GotIt, CommI
     
     IF( ParEnv % PEs<=1 ) RETURN
   
     UseL = PRESENT(LTag)
     IF(.NOT. XOR(UseL,PRESENT(Itag)) ) THEN
       CALL Fatal('CommunicateParallelSystemTag','Give either logical or integer tag!')
     END IF
     CommI = .FALSE.
     IF(.NOT. UseL) THEN
       IF(PRESENT(CommVal)) CommI = CommVal
     END IF
     
     nsize = SIZE( ParallelInfo % NodeInterface)
     IF( PRESENT(Ltag) ) THEN
       nsize = MIN(nsize, SIZE(Ltag) )
     ELSE
       nsize = MIN(nsize, SIZE(Itag) )
     END IF
     
     ALLOCATE( fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )
     
     ! Mark the neighbouring entities
     IF(ASSOCIATED( ParEnv % IsNeighbour ) ) THEN
       IsNeighbour => ParEnv % IsNeighbour
     ELSE
       ! We may want to call this even though neighbours have not been set
       ALLOCATE( IsNeighbour(ParEnv % PEs) )
       IsNeighbour = .FALSE.
       DO i=1,nsize 
         DO j=1,SIZE(ParallelInfo % Neighbourlist(i) % Neighbours)
           k = ParallelInfo % Neighbourlist(i) % Neighbours(j)
           IF ( k == ParEnv % MyPE ) CYCLE
           IsNeighbour(k+1) = .TRUE.
         END DO
       END DO
     END IF
     
     nn = 0
     ineigh = 0
     DO i=0, ParEnv % PEs-1
       k = i+1
       IF(.NOT.ParEnv % Active(k) ) CYCLE
       IF(i == ParEnv % myPE) CYCLE
       IF(.NOT. IsNeighbour(k) ) CYCLE
       nn = nn + 1
       fneigh(nn) = k
       ineigh(k) = nn
     END DO

     IF(.NOT. ASSOCIATED( ParEnv % IsNeighbour ) ) THEN
       DEALLOCATE(IsNeighbour)
     END IF
     
     ! Count the maximum number of enties to sent 
     IF( UseL ) THEN
       n = COUNT(LTag(1:nsize) .AND. ParallelInfo % NodeInterface(1:nsize))
     ELSE
       n = COUNT((ITag(1:nsize) /= 0) .AND. ParallelInfo % NodeInterface(1:nsize))
     END IF

     ! Allocate for the data to sent (s_e) and receive (r_e)
     ALLOCATE( s_e(n, nn ), r_e(n) )
     s_e = 0
     IF( CommI ) THEN
       ALLOCATE( s_i(n, nn), r_i(n) )
       s_i = 0
     END IF

     IF( CommI ) THEN
       CALL CheckBuffer( nn*6*n )
     ELSE
       CALL CheckBuffer( nn*3*n )
     END IF
       
     ii = 0
     DO i=1, nsize
       IF( UseL ) THEN
         GotIt = LTag(i) .AND. ParallelInfo % NodeInterface(i)
       ELSE
         GotIt = Itag(i) /= 0 .AND. ParallelInfo % NodeInterface(i)
       END IF
       IF(.NOT. GotIt) CYCLE
       
       DO j=1,SIZE(ParallelInfo % Neighbourlist(i) % Neighbours)
         k = ParallelInfo % Neighbourlist(i) % Neighbours(j)
         IF ( k == ParEnv % MyPE ) CYCLE
         k = k + 1
         k = ineigh(k)
         IF ( k> 0) THEN
           ii(k) = ii(k) + 1
           s_e(ii(k),k) = ParallelInfo % GlobalDOFs(i)
           IF( CommI ) THEN
             s_i(ii(k),k) = Itag(i)
           END IF
         END IF
       END DO
     END DO

     DO i=1, nn
       j = fneigh(i) 
       ! Sent size data
       CALL MPI_BSEND( ii(i),1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD,ierr )
       IF( ii(i) > 0 ) THEN
         ! Sent the global index 
         CALL MPI_BSEND( s_e(1:ii(i),i),ii(i),MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,ierr )
         IF( CommI ) THEN
           ! Sent the value of the integer tag, if requested
           CALL MPI_BSEND( s_i(1:ii(i),i),ii(i),MPI_INTEGER,j-1,112,ELMER_COMM_WORLD,ierr )
         END IF
       END IF
     END DO

     NewZeros = 0
     
     DO i=1, nn
       j = fneigh(i)
       ! Receive size of data coming from partition "j"
       CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
       IF ( n>0 ) THEN
         IF( n>SIZE(r_e)) THEN
           DEALLOCATE(r_e)
           ALLOCATE(r_e(n))
           IF(CommI) THEN
             DEALLOCATE(r_i)
             ALLOCATE(r_i(n))
           END IF
         END IF

         ! Receive the global index
         CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
         IF( CommI ) THEN
           ! Receive the value of the integer tag, if requested
           CALL MPI_RECV( r_i,n,MPI_INTEGER,j-1,112,ELMER_COMM_WORLD,status,ierr )
         END IF
         DO j=1,n
           ! Check that the entry exists in the matrix
           k = SearchNode( ParallelInfo, r_e(j), Order=ParallelInfo % Gorder )
           IF ( k>0 ) THEN
             IF( UseL ) THEN
               IF(.NOT. LTag(k)) THEN
                 LTag(k) = .TRUE.
                 NewZeros = NewZeros + 1
               END IF
             ELSE
               IF(ITag(k) == 0) THEN
                 IF( CommI ) THEN
                   ITag(k) = r_i(j)
                 ELSE
                   ITag(k) = 1
                 END IF
                 NewZeros = NewZeros + 1
               END IF
             END IF
           END IF
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e )
     IF(CommI) DEALLOCATE(s_i, r_i)

     !PRINT *,'New Zeros:',ParEnv % MyPe, NewZeros
     
  !-------------------------------------------------------------------------------
   END SUBROUTINE CommunicateParallelSystemTag
  !-------------------------------------------------------------------------------

 

 ! This subroutine fixes the global indexing of the mesh when the same mesh has been loaded to the
 ! for multiple partitions.
 !-------------------------------------------------------------------------------------------------
 SUBROUTINE SetMeshPartitionOffset(Mesh,nParMesh)
   TYPE(Mesh_t), POINTER :: Mesh  
   INTEGER :: nParMesh
   
   INTEGER :: Offset
   INTEGER :: i,n,ierr,iParExt,nParExt
   TYPE(ParallelInfo_t), POINTER :: PI

   CALL Info('SetMeshPartitionOffset','Setting offset when same mesh loaded for multiple partitions!')
   
   IF( nParMesh < 1 .OR. nParMesh >= ParEnv % PEs ) THEN
     CALL Fatal('SetMeshPartitionOffset','Invalid value of parameter nParMesh: '//TRIM(I2S(nParMesh)))
   END IF
   IF( MODULO(ParEnv % PEs, nParMesh ) /= 0 ) THEN
     CALL Fatal('SetMeshPartitionOffset','Number of partitions should be divisible with: '//TRIM(I2S(nParMesh)))
   END IF
   
   nParExt = ParEnv % PEs / nParMesh
   iParExt = ParEnv % MyPe / nParMesh

   
   PI => Mesh % ParallelInfo
   
   ! update neighbourist for partitions with an offset   
   DO i=1,Mesh % NumberOfNodes 
     IF (ASSOCIATED(PI % NeighbourList(i) % Neighbours)) THEN
       PI % NeighbourList(i) % Neighbours = &
           PI % NeighbourList(i) % Neighbours + iParExt * nParMesh
     END IF
   END DO
 
   ! Set offset for global node indexes, first find the max node index and then add the offset
   i = MAXVAL(PI % GlobalDofs )                
   CALL MPI_ALLREDUCE(i,n,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
   DO i=1,Mesh % NumberOfNodes
     PI % GlobalDofs(i) = PI % GlobalDofs(i) + iParExt * n
   END DO
   
   ! Set offset for global element indexes, first find the max element index the add the offset   
   i = MAXVAL(Mesh % Elements(:) % GElementIndex )  
   CALL MPI_ALLREDUCE(i,n,1, MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)   
   DO i=1,Mesh % NumberOfBulkElements
     Mesh % Elements(i) % GElementIndex = Mesh % Elements(i) % GElementIndex + iParExt * n
     Mesh % Elements(i) % PartIndex = Mesh % Elements(i) % PartIndex + iParExt * nParMesh
   END DO
   
 END SUBROUTINE SetMeshPartitionOffset
   
 

 SUBROUTINE InspectMesh(Mesh)
   
   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: i,j,mini,maxi
   INTEGER, POINTER :: Indexes(:)
   INTEGER, ALLOCATABLE :: ActiveCount(:)

   PRINT *,'Inspecting mesh for ranges and correctness'

   PRINT *,'No bulk elements:',Mesh % NumberOfBulkElements
   PRINT *,'No boundary elements:',Mesh % NumberOfBoundaryElements
   PRINT *,'No nodes:',Mesh % NumberOfNodes

   PRINT *,'Range:'
   PRINT *,'X:',MINVAL( Mesh % Nodes % x ), MAXVAL( Mesh % Nodes % x )
   PRINT *,'Y:',MINVAL( Mesh % Nodes % y ), MAXVAL( Mesh % Nodes % y )
   PRINT *,'Z:',MINVAL( Mesh % Nodes % z ), MAXVAL( Mesh % Nodes % z )

   ALLOCATE( ActiveCount( Mesh % NumberOfNodes ) )

   mini = HUGE(mini)
   maxi = 0
   ActiveCount = 0
   DO i=1,Mesh % NumberOfBulkElements
     Indexes => Mesh % Elements(i) % NodeIndexes
     mini = MIN(mini, MINVAL( Indexes ) )
     maxi = MAX(maxi, MAXVAL( Indexes ) )
     ActiveCount(Indexes) = ActiveCount(Indexes) + 1
   END DO
   PRINT *,'Bulk index range: ',mini,maxi
   PRINT *,'Bulk nodes:',COUNT(ActiveCount > 0 )
   PRINT *,'Bulk index count: ',MINVAL(ActiveCount),MAXVAL(ActiveCount)

   mini = HUGE(mini)
   maxi = 0
   ActiveCount = 0
   DO i=Mesh % NumberOfBulkElements+1, &
       Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
     Indexes => Mesh % Elements(i) % NodeIndexes
     mini = MIN(mini, MINVAL( Indexes ) )
     maxi = MAX(maxi, MAXVAL( Indexes ) )
     ActiveCount(Indexes) = ActiveCount(Indexes) + 1
   END DO
   PRINT *,'Boundary index range: ',mini,maxi
   PRINT *,'Boundary nodes: ',COUNT(ActiveCount > 0)
   PRINT *,'Boundary index count: ',MINVAL(ActiveCount),MAXVAL(ActiveCount)

   DEALLOCATE( ActiveCount )

   PRINT *,'Done inspecting mesh'

 END SUBROUTINE InspectMesh



!------------------------------------------------------------------------------
  SUBROUTINE SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs,NeedEdges)
!------------------------------------------------------------------------------
    INTEGER, OPTIONAL :: EdgeDOFs(:), FaceDOFs(:)
    TYPE(Mesh_t) :: Mesh
    INTEGER, OPTIONAL :: indofs(:,:)
    LOGICAL, OPTIONAL :: NeedEdges
!------------------------------------------------------------------------------
    INTEGER :: i,j,el_id
    TYPE(Element_t), POINTER :: Element, Edge, Face
    LOGICAL :: AssignEdges
!------------------------------------------------------------------------------

    CALL FindMeshEdges(Mesh)
    
    AssignEdges = .FALSE.
    IF (PRESENT(NeedEdges)) AssignEdges = NeedEdges
    
    ! Set edge and face polynomial degree and degrees of freedom for
    ! all elements
    DO i=1,Mesh % NumberOFBulkElements
       Element => Mesh % Elements(i)
       
       IF(ASSOCIATED(Element % EdgeIndexes)) THEN
         ! Iterate each edge of element
         DO j = 1,Element % TYPE % NumberOfEdges
            Edge => Mesh % Edges( Element % EdgeIndexes(j) ) 
          
            ! Set attributes of p element edges
            IF ( ASSOCIATED(Element % PDefs) ) THEN   
               ! Set edge polynomial degree and dofs
               Edge % PDefs % P = MAX( Element % PDefs % P, Edge % PDefs % P)
               Edge % BDOFs = MAX(Edge % BDOFs, Edge % PDefs % P - 1)
               Edge % PDefs % isEdge = .TRUE.
               ! Get gauss points for edge. If no dofs 2 gauss points are 
               ! still needed for integration of linear equation!
               Edge % PDefs % GaussPoints = (Edge % BDOFs+2)**Edge % TYPE % DIMENSION  

               IF (ASSOCIATED(Edge % BoundaryInfo % Left) ) THEN
                 CALL AssignLocalNumber(Edge, Edge % BoundaryInfo % Left, Mesh)
               ELSE
                 CALL AssignLocalNumber(Edge, Edge % BoundaryInfo % Right, Mesh)
               END IF
             
            ! Other element types, which need edge dofs
            ELSE IF(PRESENT(EdgeDOFs)) THEN
              Edge % BDOFs = MAX(EdgeDOFs(i), Edge % BDOFs)
            ELSE
              Edge % BDOFs = Max(1, Edge % BDOFs)
            END IF

            ! Get maximum dof for edges
            Mesh % MinEdgeDOFs = MIN(Edge % BDOFs, Mesh % MinEdgeDOFs)
            Mesh % MaxEdgeDOFs = MAX(Edge % BDOFs, Mesh % MaxEdgeDOFs)
         END DO
       END IF
       IF ( Mesh % MinEdgeDOFs > Mesh % MaxEdgeDOFs ) Mesh % MinEdgeDOFs = MEsh % MaxEdgeDOFs

       ! Iterate each face of element
       IF(.NOT. ASSOCIATED(Element % FaceIndexes)) CYCLE

       DO j=1,Element % TYPE % NumberOfFaces
          Face => Mesh % Faces( Element % FaceIndexes(j) )

          ! Set attributes of p element faces
          IF ( ASSOCIATED(Element % PDefs) ) THEN
             ! Set face polynomial degree and dofs
             Face % PDefs % P = MAX(Element % PDefs % P, Face % PDefs % P)
             ! Get number of face dofs
             Face % BDOFs = MAX( Face % BDOFs, getFaceDOFs(Element, Face % PDefs % P, j) )
             Face % PDefs % isEdge = .TRUE.
             Face % PDefs % GaussPoints = getNumberOfGaussPointsFace( Face, Mesh )
             IF (ASSOCIATED(Face % BoundaryInfo % Left) ) THEN
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Left, Mesh)
             ELSE
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Right, Mesh)
             END IF
          ELSE IF (PRESENT(FaceDOFs)) THEN
             !
             ! NOTE: This depends on what dofs have been introduced
             ! by using the construct "-quad_face b: ..." and
             ! "-tri_face b: ..."
             !
             el_id = face % TYPE % ElementCode / 100
             Face % BDOFs = MAX(FaceDOFs(i), Face % BDOFs)
             IF ( PRESENT(inDOFs) ) Face % BDOFs = MAX(Face % BDOFs, InDOFs(el_id+6,5))
          END IF
             
          ! Get maximum dof for faces
          Mesh % MinFaceDOFs = MIN(Face % BDOFs, Mesh % MinFaceDOFs)
          Mesh % MaxFaceDOFs = MAX(Face % BDOFs, Mesh % MaxFaceDOFs)
       END DO
    END DO
    IF ( Mesh % MinFaceDOFs > Mesh % MaxFaceDOFs ) Mesh % MinFaceDOFs = Mesh % MaxFaceDOFs

    ! Set local edges for boundary elements
    DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)

       ! Here set local number and copy attributes to this boundary element for left parent.
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          ! Local edges are only assigned for p elements
          IF (ASSOCIATED(Element % BoundaryInfo % Left % PDefs)) THEN
            CALL AllocatePDefinitions(Element)
            Element % PDefs % isEdge = .TRUE.
            CALL AssignLocalNumber(Element, Element % BoundaryInfo % Left, Mesh)
            ! CYCLE
          END IF
       END IF

       ! Here set local number and copy attributes to this boundary element for right parent
       IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          ! Local edges are only assigned for p elements
          IF (ASSOCIATED(Element % BoundaryInfo % Right % PDefs)) THEN
             CALL AllocatePDefinitions(Element)
             Element % PDefs % isEdge = .TRUE.
             CALL AssignLocalNumber(Element, Element % BoundaryInfo % Right, Mesh)
          END IF
       END IF

       IF (AssignEdges) THEN
         IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
           CALL AssignLocalNumber(Element,Element % BoundaryInfo % Left, Mesh, NoPE=.TRUE.)
         END IF
         IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
           CALL AssignLocalNumber(Element,Element % BoundaryInfo % Right, Mesh, NoPE=.TRUE.)
         END IF
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetMeshEdgeFaceDofs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE SetMeshMaxDOFs(Mesh)
!------------------------------------------------------------------------------
   TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
   TYPE(Element_t), POINTER :: Element
   INTEGER :: i,j,n

   DO i=1,Mesh % NumberOfBulkElements
     Element => Mesh % Elements(i)

     ! Set gauss points for each p element
     IF ( ASSOCIATED(Element % PDefs) ) THEN
       Element % PDefs % GaussPoints = getNumberOfGaussPoints( Element, Mesh )
     END IF

     Mesh % MaxBDOFs = MAX( Element % BDOFs, Mesh % MaxBDOFs )
     Mesh % MaxNDOFs = MAX(Element % NDOFs / Element % TYPE % NumberOfNodes, &
         Mesh % MaxNDOFs)
   END DO

   DO i=1,Mesh % NumberOFBulkElements
     Element => Mesh % Elements(i)

     ! Set max element dofs here (because element size may have changed
     ! when edges and faces have been set). This is the absolute worst case.
     ! Element which has MaxElementDOFs may not even be present as a 
     ! real element
     Mesh % MaxElementDOFs = MAX( Mesh % MaxElementDOFs, &
          Element % TYPE % NumberOfNodes * Mesh % MaxNDOFs + &
          Element % TYPE % NumberOfEdges * Mesh % MaxEdgeDOFs + &
          Element % TYPE % NumberOfFaces * Mesh % MaxFaceDOFs + &
          Element % BDOFs, &
          Element % DGDOFs )

     IF ( Element % BDOFs > 0 ) THEN
       ALLOCATE( Element % BubbleIndexes(Element % BDOFs) )
       DO j=1,Element % BDOFs
         Element % BubbleIndexes(j) = Mesh % MaxBDOFs*(i-1)+j
       END DO
     END IF
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE SetMeshMaxDOFs
!------------------------------------------------------------------------------
 
 SUBROUTINE ReadTargetNames(Model,Filename)
     CHARACTER(LEN=*) :: FileName
     TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
   INTEGER, PARAMETER :: FileUnit = 10
   INTEGER, PARAMETER :: A=ICHAR('A'),Z=ICHAR('Z'),U2L=ICHAR('a')-ICHAR('A')
   INTEGER :: i,j,k,iostat,i1,i2,i3,n
   INTEGER :: ivals(256)
   CHARACTER(LEN=1024) :: str, name0, name1
   TYPE(ValueList_t), POINTER :: Vlist
   LOGICAL :: Found, AlreadySet

   OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT=iostat )
   IF( iostat /= 0 ) THEN
     CALL Fatal('ReadTargetNames','Requested the use of entity names but this file does not exits: '//TRIM(FileName))
   END IF
   
   CALL Info('ReadTargetNames','Reading names info from file: '//TRIM(FileName))

   DO WHILE( .TRUE. ) 
     READ(FileUnit,'(A)',IOSTAT=iostat) str
     IF( iostat /= 0 ) EXIT
     i = INDEX( str,'$')     
     j = INDEX( str,'=')
     IF( i == 0 .OR. j == 0 ) CYCLE

     i = i + 1
     DO WHILE(i<=LEN_TRIM(str) .AND. str(i:i)==' ')
       i = i + 1
     END DO     
     
     i1 = i
     i2 = j-1
     i3 = j+1

     ! Move to lowercase since the "name" in sif file is also
     ! always in lowercase. 
     DO i=i1,i2
       j = i+1-i1
       k = ICHAR(str(i:i))
       IF ( k >= A .AND. k<= Z ) THEN
         name0(j:j) = CHAR(k+U2L)
       ELSE
         name0(j:j) = str(i:i)
       END IF
     END DO

     n = str2ints( str(i3:),ivals )
     IF( n == 0 ) THEN
       CALL Fatal('ReadTargetNames','Could not find arguments for: '//str(i1:i2))
     END IF

     AlreadySet = .FALSE.

     DO i=1,Model % NumberOfBCs
       Vlist => Model % BCs(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches BC '//TRIM(I2S(i))
         IF( AlreadySet ) THEN
           CALL Fatal('ReadTargetNames','Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Boundaries') ) THEN
           CALL Info('ReadTargetNames','> Target Boundaries < already defined for BC '&
               //TRIM(I2S(i)))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Boundaries',n,ivals(1:n))
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO

     DO i=1,Model % NumberOfBodies
       Vlist => Model % Bodies(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches body '//TRIM(I2S(i))
         IF( AlreadySet ) THEN
           CALL Fatal('ReadTargetNames','Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Bodies') ) THEN
           CALL Info('ReadTargetNames','> Target Bodies < already defined for Body '&
               //TRIM(I2S(i)))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Bodies',n,ivals(1:n))
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO
     
     IF(.NOT. AlreadySet ) THEN
       CALL Warn('ReadTargetNames','Could not map name to Body nor BC: '//name0(1:i2-i1+1) )
     END IF

   END DO

   CLOSE(FileUnit)
   
 END SUBROUTINE ReadTargetNames


!------------------------------------------------------------------------------
!> This subroutine reads elementwise input data from the file mesh.elements.data 
!> and inserts the data into the structured data variable 
!> Mesh % Elements(element_id) % PropertyData. The contents of the file should
!> be arranged as
!> 
!> element: element_id_1
!> data_set_name_1: a_1 a_2 ... a_n
!> data_set_name_2: b_1 b_2 ... b_m
!> data_set_name_3: ...
!> end
!> element: ...
!> ...
!> end
!------------------------------------------------------------------------------
  SUBROUTINE ReadElementPropertyFile(FileName,Mesh)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: FileName
     TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MAXLEN=1024
    CHARACTER(LEN=:), ALLOCATABLE :: str
    INTEGER :: i,j,n
    INTEGER, PARAMETER :: FileUnit = 10
    REAL(KIND=dp) :: x
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementData_t), POINTER :: PD,PD1
!------------------------------------------------------------------------------
    ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)

    OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', ERR=10 )

    DO WHILE( ReadAndTrim(FileUnit,str) )
      READ( str(9:),*) i
      IF ( i < 0 .OR. i > Mesh % NumberOFBulkElements ) THEN
        CALL Fatal( 'ReadElementPropertyFile', 'Element id out of range.' )
      END IF

      IF ( SEQL( str, 'element:') ) THEN
        Element => Mesh % Elements(i)
        PD => Element % PropertyData

        DO WHILE(ReadAndTrim(FileUnit,str))
          IF ( str == 'end' ) EXIT

          i = INDEX(str, ':')
          IF ( i<=0 ) CYCLE

          IF ( .NOT.ASSOCIATED(PD)  ) THEN
            ALLOCATE( Element % PropertyData )
            PD => Element % PropertyData
            PD % Name = TRIM(str(1:i-1))
          ELSE
            DO WHILE(ASSOCIATED(PD))
              IF ( PD % Name==TRIM(str(1:i-1)) ) EXIT
              PD1 => PD
              PD => PD % Next
            END DO
            
            IF (.NOT. ASSOCIATED(PD) ) THEN
              ALLOCATE(PD1 % Next)
              PD => PD1 % Next
              PD % Name = TRIM(str(1:i-1))
            END IF
          END IF

          j = i+1
          n = 0
          DO WHILE(j<=LEN_TRIM(str))
            READ( str(j:), *, END=20,ERR=20 ) x
            n = n + 1
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
              j = j + 1
            END DO
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
              j = j + 1
            END DO
          END DO
20        CONTINUE
          IF ( n>0 ) THEN
            ALLOCATE(PD % Values(n))
            j = i+1
            n = 1
            DO WHILE(j<=LEN_TRIM(str))
              READ( str(j:), *, END=30,ERR=30 ) PD % Values(n)
              n = n + 1
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
                j = j + 1
              END DO
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
                j = j + 1
              END DO
            END DO
30          CONTINUE
          END IF
        END DO
      END IF
    END DO

    CLOSE(FileUnit)

10  CONTINUE

!------------------------------------------------------------------------------
  END SUBROUTINE ReadElementPropertyFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE MeshStabParams( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    INTEGER :: i,n, istat
    LOGICAL :: stat, Stabilize, UseLongEdge
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

    CALL Info('MeshStabParams','Computing stabilization parameters',Level=7)
    CALL ResetTimer('MeshStabParams')

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('MeshStabParams','Mesh not associated')
    END IF
    
    IF ( Mesh % NumberOfNodes <= 0 ) RETURN

    Stabilize = .FALSE.
    
    DO i=1,CurrentModel % NumberOfSolvers
      Solver => CurrentModel % Solvers(i)
      IF ( ASSOCIATED( Mesh, Solver % Mesh ) ) THEN
        Stabilize = Stabilize .OR. &
            ListGetLogical( Solver % Values, 'Stabilize', Stat )
        Stabilize = Stabilize .OR. &
            ListGetString( Solver % Values,  &
            'Stabilization Method', Stat )=='vms'
        Stabilize = Stabilize .OR. &
            ListGetString( Solver % Values,  &
            'Stabilization Method', Stat )=='stabilized'
      END IF
    END DO

    Mesh % Stabilize = Stabilize 
    
    IF( ListGetLogical(CurrentModel % Simulation, &
        "Skip Mesh Stabilization",Stat) ) RETURN
    
    !IF( .NOT. Stabilize ) THEN
    !  CALL Info('MeshStabParams','No need to compute stabilization parameters',Level=10)      
    !  RETURN      
    !END IF
    
    CALL AllocateVector( Nodes % x, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % y, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % z, Mesh % MaxElementNodes )

    UseLongEdge = ListGetLogical(CurrentModel % Simulation, &
         "Stabilization Use Longest Element Edge",Stat)

    DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       n = Element % TYPE % NumberOfNodes
       Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
       Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
       Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
       IF ( Mesh % Stabilize ) THEN
          CALL StabParam( Element, Nodes,n, &
              Element % StabilizationMK, Element % hK, UseLongEdge=UseLongEdge)
       ELSE
          Element % hK = ElementDiameter( Element, Nodes, UseLongEdge=UseLongEdge)
       END IF
    END DO
 
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

    CALL CheckTimer('MeshStabParams',Level=7,Delete=.TRUE.)
!----------------------------------------------------------------------------
  END SUBROUTINE MeshStabParams
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The quadratic mesh should be such that the center nodes lie roughly between
!> the corner nodes. This routine checks that this is actually the case.
!> The intended use for the routine is different kind of mesh related debugging.
!------------------------------------------------------------------------------
  SUBROUTINE InspectQuadraticMesh( Mesh, EnforceToCenter ) 
    
    TYPE(Mesh_t), TARGET :: Mesh
    LOGICAL, OPTIONAL :: EnforceToCenter

    LOGICAL :: Enforce
    INTEGER :: i,n,k,k1,k2,k3,ElemCode,ElemFamily,ElemDegree,ErrCount,TotCount
    REAL(KIND=dp) :: Center(3),Ref(3),Dist,Length
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: CenterMap(:,:)
    INTEGER, TARGET  :: TriangleCenterMap(3,3), QuadCenterMap(4,3), &
        TetraCenterMap(6,3), BrickCenterMap(12,3), WedgeCenterMap(9,3), PyramidCenterMap(8,3) 
    
    CALL Info('InspectQuadraticMesh','Inspecting quadratic mesh for outliers')
    CALL Info('InspectQuadraticMesh','Number of nodes: '//TRIM(I2S(Mesh % NumberOfNodes)),Level=8)
    CALL Info('InspectQuadraticMesh','Number of bulk elements: '&
        //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=8)
    CALL Info('InspectQuadraticMesh','Number of boundary elements: '&
        //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=8)


    IF( PRESENT( EnforceToCenter ) ) THEN
      Enforce = EnforceToCenter
    ELSE
      Enforce = .FALSE.
    END IF

    TriangleCenterMap(1,:) = [ 1, 2, 4]
    TriangleCenterMap(2,:) = [ 2, 3, 5]
    TriangleCenterMap(3,:) = [ 3, 1, 6]
    
    QuadCenterMap(1,:) = [ 1, 2, 5]
    QuadCenterMap(2,:) = [ 2, 3, 6]
    QuadCenterMap(3,:) = [ 3, 4, 7]
    QuadCenterMap(4,:) = [ 4, 1, 8]
    
    TetraCenterMap(1,:) = [ 1, 2, 5]
    TetraCenterMap(2,:) = [ 2, 3, 6]
    TetraCenterMap(3,:) = [ 3, 1, 7]
    TetraCenterMap(4,:) = [ 1, 4, 8]
    TetraCenterMap(5,:) = [ 2, 4, 9]
    TetraCenterMap(6,:) = [ 3, 4, 10]

    BrickCenterMap(1,:) = [ 1, 2,  9 ]
    BrickCenterMap(2,:) = [ 2, 3,  10 ]
    BrickCenterMap(3,:) = [ 3, 4,  11 ]
    BrickCenterMap(4,:) = [ 4, 1,  12 ]
    BrickCenterMap(5,:) = [ 1, 5,  13 ]
    BrickCenterMap(6,:) = [ 2, 6,  14 ]
    BrickCenterMap(7,:) = [ 3, 7,  15 ]
    BrickCenterMap(8,:) = [ 4, 8,  16 ]
    BrickCenterMap(9,:) = [ 5, 6,  17 ]
    BrickCenterMap(10,:) = [ 6, 7, 18 ]
    BrickCenterMap(11,:) = [ 7, 8, 19 ]
    BrickCenterMap(12,:) = [ 8, 5, 20 ]
    
    WedgeCenterMap(1,:) = [ 1, 2, 7 ]
    WedgeCenterMap(2,:) = [ 2, 3, 8 ]
    WedgeCenterMap(3,:) = [ 3, 1, 9 ]
    WedgeCenterMap(4,:) = [ 4, 5, 10 ]
    WedgeCenterMap(5,:) = [ 5, 6, 11 ]
    WedgeCenterMap(6,:) = [ 6, 4, 12 ]
    WedgeCenterMap(7,:) = [ 1, 4, 13 ]
    WedgeCenterMap(8,:) = [ 2, 5, 14 ]
    WedgeCenterMap(9,:) = [ 3, 6, 15 ]
    
    PyramidCenterMap(1,:) = [ 1,2,6 ]
    PyramidCenterMap(2,:) = [ 2,3,7 ]
    PyramidCenterMap(3,:) = [ 3,4,8 ]
    PyramidCenterMap(4,:) = [ 4,1,9 ]
    PyramidCenterMap(5,:) = [ 1,5,10 ]
    PyramidCenterMap(6,:) = [ 2,5,11 ]
    PyramidCenterMap(7,:) = [ 3,5,12 ]
    PyramidCenterMap(8,:) = [ 4,5,13 ]
    
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z
    
    !   Loop over elements:
    !   -------------------
    ErrCount = 0
    TotCount = 0

    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(i)

      ElemCode = Element % TYPE % ElementCode 
      ElemFamily = ElemCode / 100
      ElemDegree = Element % TYPE % BasisFunctionDegree
      
      ! Only check quadratic elements!
      IF( ElemDegree /= 2 ) CYCLE
      
      SELECT CASE( ElemFamily ) 

      CASE(3)
        n = 3
        CenterMap => TriangleCenterMap
        
      CASE(4)
        n = 4
        CenterMap => QuadCenterMap
        
      CASE(5)
        n = 6
        CenterMap => TetraCenterMap
        
      CASE(6)
        n = 8
        CenterMap => PyramidCenterMap
        
      CASE(7)
        n = 9
        CenterMap => WedgeCenterMap
        
      CASE(8)
        n = 12
        CenterMap => BrickCenterMap
        
      CASE DEFAULT
        CALL Fatal('InspectQuadraticMesh','Element type '//TRIM(I2S(ElemCode))//' not implemented!')

      END SELECT
      
      !      Loop over every edge of every element:
      !      --------------------------------------
       DO k=1,n
         k1 = Element % NodeIndexes( CenterMap(k,1) )
         k2 = Element % NodeIndexes( CenterMap(k,2) )
         k3 = Element % NodeIndexes( CenterMap(k,3) )
         
         Center(1) = ( x(k1) + x(k2) ) / 2.0_dp
         Center(2) = ( y(k1) + y(k2) ) / 2.0_dp
         Center(3) = ( z(k1) + z(k2) ) / 2.0_dp

         Ref(1) = x(k3)
         Ref(2) = y(k3) 
         Ref(3) = z(k3)

         Length = SQRT( (x(k1) - x(k2))**2.0 + (y(k1) - y(k2))**2.0 + (z(k1) - z(k2))**2.0 )
         Dist = SQRT( SUM( (Center - Ref)**2.0 ) )

         TotCount = TotCount + 1
         IF( Dist > 0.01 * Length ) THEN
           ErrCount = ErrCount + 1
           PRINT *,'Center Displacement:',i,ElemCode,n,k,Dist/Length
         END IF

         IF( Enforce ) THEN
           x(k3) = Center(1)
           y(k3) = Center(2)
           z(k3) = Center(3)
         END IF

       END DO
     END DO
         
     IF( TotCount > 0 ) THEN
       CALL Info('InspectQuadraticMesh','Number of outlier nodes is '&
           //TRIM(I2S(ErrCount))//' out of '//TRIM(I2S(TotCount)),Level=6)
     ELSE
       CALL Info('InspectQuadraticMesh','No quadratic elements to inspect',Level=8)
     END IF

  END SUBROUTINE InspectQuadraticMesh
 

  
  SUBROUTINE FollowCurvedBoundary(Model, Mesh, SetP )
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh 
    LOGICAL, OPTIONAL :: SetP

    LOGICAL :: Found
    REAL(KIND=dp) :: FitParams(7)
    INTEGER :: Mode, bc_ind, dim
    TYPE(ValueList_t), POINTER :: Vlist
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshName
    TYPE(Mesh_t), POINTER :: BMesh

    IF(.NOT. ListCheckPrefixAnyBC( Model,'Follow') .OR. &
        ListCheckPrefix( Model % Simulation,'Follow') ) RETURN

    dim = Mesh % MeshDim
    
    DO bc_ind = 0, Model % NumberOfBCs
      IF( bc_ind == 0 ) THEN
        Vlist => Model % Simulation
      ELSE
        Vlist => Model % BCs(bc_ind) % Values
      END IF      

      Mode = 0
      IF( ListGetLogical(Vlist,'Follow Circle Boundary', Found ) ) THEN
        CALL CylinderFit(Mesh, Vlist, bc_ind, 2, FitParams ) 
        Mode = 1        
      ELSE IF( ListGetLogical(Vlist,'Follow Cylinder Boundary', Found ) ) THEN
        CALL CylinderFit(Mesh, Vlist, bc_ind, dim, FitParams) 
        Mode = 2        
      ELSE IF( ListGetLogical(Vlist,'Follow Sphere Boundary', Found ) ) THEN
        CALL SphereFit(Mesh, Vlist, bc_ind, FitParams ) 
        Mode = 3        
      ELSE IF( ListGetLogical(Vlist,'Follow Function Boundary', Found ) ) THEN
        IF(.NOT. ListCheckPresent(Vlist,'Surface Function') ) THEN
          CALL Fatal('FollowCurvedBoundary','We need "Surface Function" to follow!')
        END IF
        Mode = 4        
      ELSE IF( ListCheckPresent(Vlist,'Follow Mesh') ) THEN
        MeshName = ListGetString( Vlist,'Follow Mesh' )
        Found = .FALSE.
        BMesh => Model % Meshes
        DO WHILE( ASSOCIATED( BMesh ) )          
          IF(ASSOCIATED(Bmesh,Mesh) ) BMesh => BMesh % Next
          IF(TRIM(MeshName) == TRIM(BMesh % Name) ) THEN
            Found = .TRUE.
            EXIT
          END IF
        END DO
        IF(.NOT. Found ) THEN
          CALL Fatal('FollowCurvedBoundary','Could not find mesh to follow: '//TRIM(MeshName))
        END IF
        Mode = 5
      END IF
      
      IF(Mode > 0 ) THEN
        IF(bc_ind == 0 ) THEN
          CALL Info('FollowCurvedBoundary','Setting whole mesh '//&
              ' to follow curved boundary in mode '//TRIM(I2S(Mode)),Level=7)
        ELSE
          CALL Info('FollowCurvedBoundary','Setting BC '//TRIM(I2S(bc_ind))//&
              ' to follow curved boundary in mode '//TRIM(I2S(Mode)),Level=7)
        END IF
        CALL SetCurvedBoundary()
      END IF
    END DO

    
  CONTAINS

          
!------------------------------------------------------------------------------
    SUBROUTINE SetCurvedBoundary()
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: R, rat, f, gradf(3)
      REAL(KIND=dp) :: Nrm(3), Tngt1(3), Tngt2(3), Orig(3), Coord(3), NtCoord(3)
      INTEGER :: i,j,k,l,t,n,t1,t2
      LOGICAL, POINTER :: DoneNode(:)
      TYPE(Element_t), POINTER :: Element
      LOGICAL :: Parallel 
      TYPE(ParallelInfo_t), POINTER :: ParallelInfo
      INTEGER, TARGET :: Hdim(1)
      INTEGER, POINTER :: pHdim(:)
      
      
      IF( Mode == 1 ) THEN  ! circle
        Orig(1:2) = FitParams(1:2)
        Orig(3) = 0.0_dp
        R = FitParams(3)
        IF( InfoActive(25) .AND. ParEnv % MyPe == 0) PRINT *,'Circle Params:',FitParams(1:3)                        
      ELSE IF( Mode == 2 ) THEN  ! cylinder 
        Orig(1:3) = FitParams(1:3)
        Nrm(1:3) = FitParams(4:6)        
        R = FitParams(7)
        IF( InfoActive(25) .AND. ParEnv % MyPe == 0) PRINT *,'Cylinder Params:',FitParams(1:7)        
        CALL TangentDirections(Nrm, Tngt1, Tngt2 ) 
      ELSE IF( Mode == 3 ) THEN ! sphere
        Orig(1:3) = FitParams(1:3)
        Nrm = 0.0_dp
        R = FitParams(4)
        IF( InfoActive(25) .AND. ParEnv % MyPe == 0) PRINT *,'Sphere Params:',FitParams(1:4)                                
      ELSE IF( Mode == 4 ) THEN
        Orig = 0.0_dp        
      END IF
      
      Parallel = ( ParEnv % PEs > 1 .AND. .NOT. Mesh % SingleMesh )
      IF( bc_ind == 0 ) THEN
        t1 = 1
        t2 = Mesh % NumberOfBulkElements
      ELSE
        t1 = Mesh % NumberOfBulkElements + 1
        t2 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      END IF
      
      IF(.NOT. SetP) THEN
        
        IF( bc_ind > 0 ) THEN
          ALLOCATE( DoneNode(Mesh % NumberOfNodes))
          DoneNode = .FALSE.
        
          DO t=t1, t2
            Element => Mesh % Elements(t)
            IF ( Element % BoundaryInfo % Constraint &
                /= Model % BCs(bc_ind) % Tag ) CYCLE
            n = Element % TYPE % NumberOfNodes          
            DoneNode(Element % NodeIndexes(1:n)) = .TRUE.
          END DO
          
          IF( Parallel ) THEN
            ParallelInfo => Mesh % ParallelInfo 
            CALL CommunicateParallelSystemTag(ParallelInfo,Ltag = DoneNode)
          END IF
        END IF

        IF( Mode == 5 ) THEN          
            Hdim(1) = 3
            pHdim => Hdim
            IF( bc_ind > 0 ) THEN          
              CALL InterpolateVartoVarReduced( BMesh, Mesh,'Coordinate 3', pHdim,NewNodeMask = DoneNode )
            ELSE
              CALL InterpolateVartoVarReduced( BMesh, Mesh,'Coordinate 3', pHdim )
          END IF                    
          RETURN
        END IF

        
        DO j=1, Mesh % NumberOfNodes
          IF( Bc_ind > 0 ) THEN
            IF( .NOT. DoneNode(j) ) CYCLE
          END IF
            
          Coord(1) = Mesh % Nodes % x(j) - Orig(1)           
          Coord(2) = Mesh % Nodes % y(j) - Orig(2)
          Coord(3) = Mesh % Nodes % z(j) - Orig(3)
          
          SELECT CASE( Mode )
          CASE( 1 ) ! circle 
            rat = R / SQRT(SUM(Coord(1:2)**2))
            Coord(1:2) = rat*Coord(1:2)
          CASE( 2 ) ! cylinder
            NtCoord(1) = SUM(Nrm*Coord)
            NtCoord(2) = SUM(Tngt1*Coord)
            NtCoord(3) = SUM(Tngt2*Coord)
            rat = R / SQRT(SUM(NtCoord(2:3)**2))
            NtCoord(2:3) = rat*NtCoord(2:3)
            Coord = NtCoord(1)*Nrm + NtCoord(2)*Tngt1 + NtCoord(3)*Tngt2
          CASE( 3 ) ! sphere 
            rat = R / SQRT(SUM(Coord(1:3)**2))
            Coord(1:3) = rat*Coord(1:3)
          CASE( 4 ) ! analytical function
            ! For now we fix Newton's iteration to three...
            DO i=1,3
              f = ListGetFunVec( Vlist,'Surface Function', Coord(1:dim), dim, DfDx=gradf(1:dim) )
              Coord(1:dim) = Coord(1:dim) - f*gradf(1:dim)/(SUM(gradf(1:dim)**2))
            END DO
          END SELECT
          
          Mesh % Nodes % x(j) = Coord(1) + Orig(1)
          Mesh % Nodes % y(j) = Coord(2) + Orig(2)
          Mesh % Nodes % z(j) = Coord(3) + Orig(3)
        END DO
        IF( Bc_ind > 0 ) THEN
          DEALLOCATE(DoneNode)
        END IF
      END IF
        
      IF( SetP ) THEN
        DO t=t1, t2
          Element => Mesh % Elements(t)
          IF( bc_ind > 0 ) THEN
            IF ( Element % BoundaryInfo % Constraint &
                /= Model % BCs(bc_ind) % Tag ) CYCLE
          END IF
          n = Element % TYPE % NumberOfNodes
          
          BLOCK 
            REAL(KIND=dp) :: Weight
            REAL(KIND=dp) :: Basis(50),DetJ
            REAL(KIND=dp) :: MASS(50,50), FORCE(3,50), x(50), Coord0(3)
            LOGICAL :: Stat, Erroneous
            INTEGER :: nd,i,t,p,q
            INTEGER, TARGET :: Indexes(50)
            INTEGER :: pivot(50)
            INTEGER, POINTER :: pIndexes(:)
            TYPE(GaussIntegrationPoints_t) :: IP
            TYPE(Nodes_t), SAVE :: Nodes

            pIndexes => Indexes 
            Nd = mGetElementDOFs( pIndexes, Element, CurrentModel % Solver )          
                        
            ! Only if we have really p-elements is there a need to consider the curved shape
            IF(Nd == n ) CYCLE

            CALL CopyElementNodesFromMesh( Nodes, Mesh, n, pIndexes)

            MASS = 0._dp
            FORCE = 0._dp

            IP = GaussPoints( Element )
            
            DO t=1,IP % n
              stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis )
              Weight = IP % s(t) * DetJ

              ! Current nodal value at integration point does not consider p-dofs
              Coord(1) = SUM( Nodes % x(1:n) * Basis(1:n) )
              Coord(2) = SUM( Nodes % y(1:n) * Basis(1:n) )
              Coord(3) = SUM( Nodes % z(1:n) * Basis(1:n) )
              Coord0 = Coord

              Coord = Coord - Orig
              SELECT CASE( Mode )
              CASE( 1 ) 
                rat = R / SQRT(SUM(Coord(1:2)**2))
                Coord(1:2) = rat * Coord(1:2)
              CASE( 2 )
                ! Local coordinates in nt-system
                NtCoord(1) = SUM(Nrm*Coord)
                NtCoord(2) = SUM(Tngt1*Coord)
                NtCoord(3) = SUM(Tngt2*Coord)
                ! Ratio between current and desired radius
                rat = R / SQRT(SUM(NtCoord(2:3)**2))
                NtCoord(2:3) = rat * NtCoord(2:3)
                Coord = NtCoord(1)*Nrm + NtCoord(2)*Tngt1 + NtCoord(3)*Tngt2
              CASE( 3 ) 
                rat = R / SQRT(SUM(Coord(1:3)**2))
                Coord(1:3) = rat * Coord(1:3)
              CASE( 4 ) 
                DO i=1,3
                  f = ListGetFunVec( Vlist,'Surface Function', Coord(1:dim), dim, DfDx=gradf(1:dim) )
                  Coord(1:dim) = Coord(1:dim) - f*gradf(1:dim)/(SUM(gradf(1:dim)**2))            
                END DO
              END SELECT

              Coord = Coord + Orig
              ! Solve for desired coordinate displacement rather than absolute coordinate value
              Coord = Coord - Coord0
                
              ! Create equation involving mass matrix that solves for the coordinates at the p-dofs
              DO q=1,nd
                MASS(1:nd,q) = MASS(1:nd,q) + Weight * Basis(1:nd) * Basis(q) 
              END DO

              DO i=1,dim
                FORCE(i,1:nd) = FORCE(i,1:nd) + Weight * Basis(1:nd) * Coord(i) 
              END DO
            END DO

            ! Set Dirichlet conditions for the nodal coordinate displacements
            DO i=1,n
              MASS(i,1:nd) = 0.0_dp
              MASS(i,i) = 1.0_dp
              FORCE(:,i) = 0.0_dp
            END DO
            
            CALL LUdecomp(MASS,nd,pivot,Erroneous)
            IF (Erroneous) CALL Fatal('SetCurvedBoundary', 'LU-decomposition fails')
            
            DO i=1,dim          
              x(1:nd) = FORCE(i,1:nd)
              CALL LUSolve(nd,MASS,x,pivot)
              
              SELECT CASE(i)
              CASE(1)
                Mesh % Nodes % x(Indexes(n+1:nd)) = x(n+1:nd) 
              CASE(2)
                Mesh % Nodes % y(Indexes(n+1:nd)) = x(n+1:nd) 
              CASE(3)
                Mesh % Nodes % z(Indexes(n+1:nd)) = x(n+1:nd) 
              END SELECT
            END DO
            
          END BLOCK
        END DO
      END IF
        
    END SUBROUTINE SetCurvedBoundary
!------------------------------------------------------------------------------
  END SUBROUTINE FollowCurvedBoundary

  
  
  !------------------------------------------------------------------------------------------------
  !> Finds nodes for which CandNodes are True such that their mutual distance is somehow
  !> maximized. We first find lower left corner, then the node that is furtherst apart from it,
  !> and continue as long as there are nodes to find. Typically we would be content with two nodes
  !> on a line, three nodes on a plane, and four nodes on a volume.
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE FindExtremumNodes(Mesh,CandNodes,NoExt,Inds) 
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, ALLOCATABLE :: CandNodes(:)
    INTEGER :: NoExt
    INTEGER, POINTER :: Inds(:)

    REAL(KIND=dp) :: Coord(3),dCoord(3),dist,MinDist,MaxDist
    REAL(KIND=dp), ALLOCATABLE :: SetCoord(:,:)
    INTEGER :: i,j,k
    
    ALLOCATE( SetCoord(NoExt,3) )
    SetCoord = 0.0_dp
    Inds = 0
    
    ! First find the lower left corner
    MinDist = HUGE(MinDist) 
    DO i=1, Mesh % NumberOfNodes
      IF(.NOT. CandNodes(i) ) CYCLE
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
      Dist = SUM( Coord )
      IF( Dist < MinDist ) THEN
        Inds(1) = i
        MinDist = Dist
        SetCoord(1,:) = Coord
      END IF
    END DO
    
    ! Find more points such that their minimum distance to the previous point(s)
    ! is maximized.
    DO j=2,NoExt
      ! The maximum minimum distance of any node from the previously defined nodes
      MaxDist = 0.0_dp
      DO i=1, Mesh % NumberOfNodes
        IF(.NOT. CandNodes(i) ) CYCLE
        Coord(1) = Mesh % Nodes % x(i)
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)
        
        ! Minimum distance from the previously defined nodes
        MinDist = HUGE(MinDist)
        DO k=1,j-1
          dCoord = SetCoord(k,:) - Coord
          Dist = SUM( dCoord**2 )          
          MinDist = MIN( Dist, MinDist )
        END DO
        
        ! If the minimum distance is greater than in any other node, choose this
        IF( MaxDist < MinDist ) THEN
          MaxDist = MinDist 
          Inds(j) = i
          SetCoord(j,:) = Coord
        END IF
      END DO
    END DO

    IF( InfoActive(20) ) THEN
      PRINT *,'Extremum Inds:',Inds
      DO i=1,NoExt
        PRINT *,'Node:',Inds(i),SetCoord(i,:)
      END DO
    END IF
      
  END SUBROUTINE FindExtremumNodes
    

  

  ! Save projector, mainly a utility for debugging purposes
  !--------------------------------------------------------
  SUBROUTINE SaveProjector(Projector,SaveRowSum,Prefix,InvPerm,Parallel)
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: SaveRowSum 
    CHARACTER(LEN=*) :: Prefix
    INTEGER, POINTER, OPTIONAL :: InvPerm(:)
    LOGICAL, OPTIONAL :: Parallel

    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    INTEGER :: i,j,ii,jj
    REAL(KIND=dp) :: rowsum, dia, val
    INTEGER, POINTER :: IntInvPerm(:)
    LOGICAL :: GlobalInds
    INTEGER, POINTER :: GlobalDofs(:)
    CHARACTER(*), PARAMETER :: Caller = "SaveProjector"
    
    IF(.NOT.ASSOCIATED(Projector)) RETURN
    
    IF( PRESENT( InvPerm ) ) THEN
      IntInvPerm => InvPerm 
    ELSE
      IntInvPerm => Projector % InvPerm
    END IF

    GlobalInds = .FALSE.
    IF(ParEnv % PEs == 1 ) THEN
      FileName = TRIM(Prefix)//'.dat'
    ELSE
      FileName = TRIM(Prefix)//'_part'//&
          TRIM(I2S(ParEnv % MyPe))//'.dat'
      IF( PRESENT( Parallel ) ) GlobalInds = Parallel
    END IF

    IF( GlobalInds ) THEN
      NULLIFY( GlobalDofs ) 
      IF( ASSOCIATED( CurrentModel % Solver % Matrix ) ) THEN
        GlobalDofs => CurrentModel % Solver % Matrix % ParallelInfo % GlobalDofs
      END IF
      IF(.NOT. ASSOCIATED( GlobalDofs ) ) THEN
        CALL Info(Caller,'Cannot find GlobalDofs for Solver matrix')
        GlobalDofs => CurrentModel % Mesh % ParallelInfo % GlobalDofs
      END IF
    END IF
          
    OPEN(1,FILE=FileName,STATUS='Unknown')    
    DO i=1,projector % numberofrows
      IF( ASSOCIATED( IntInvPerm ) ) THEN
        ii = intinvperm(i)        
        IF( ii == 0) THEN
          PRINT *,'Projector InvPerm is zero:',ParEnv % MyPe, i, ii
          CYCLE
        END IF
      ELSE
        ii = i
      END IF
      IF( GlobalInds ) THEN
        IF( ii > SIZE( GlobalDofs ) ) THEN
          PRINT *,'ParEnv % MyPe, Projecor invperm is larger than globaldofs',&
              ii, SIZE( GlobalDofs ), i, Projector % NumberOfRows
          CYCLE
        END IF
        ii = GlobalDofs(ii)
      END IF
      IF( ii == 0) THEN
        PRINT *,'Projector global InvPerm is zero:',ParEnv % MyPe, i, ii
        CYCLE
      END IF
      DO j=projector % rows(i), projector % rows(i+1)-1
        jj = projector % cols(j)
        IF( jj == 0) THEN
          PRINT *,'Projector col is zero:',ParEnv % MyPe, i, ii, j, jj
          CYCLE
        END IF       
        val = projector % values(j)
        IF( GlobalInds ) THEN
          IF( jj > SIZE( GlobalDofs ) ) THEN
            PRINT *,'Projecor invperm is larger than globaldofs',&
                jj, SIZE( GlobalDofs )
            CYCLE
          END IF
          jj = GlobalDofs(jj)
          IF( jj == 0) THEN
            PRINT *,'Projector global col is zero:',ParEnv % MyPe, i, ii, j, jj
            CYCLE
          END IF
          WRITE(1,*) ii,jj,ParEnv % MyPe, val
        ELSE
          WRITE(1,*) ii,jj,val
        END IF
      END DO
    END DO
    CLOSE(1)     

    IF( SaveRowSum ) THEN
      IF(ParEnv % PEs == 1 ) THEN
        FileName = TRIM(Prefix)//'_rsum.dat'
      ELSE
        FileName = TRIM(Prefix)//'_rsum_part'//&
            TRIM(I2S(ParEnv % MyPe))//'.dat'
      END IF
      
      OPEN(1,FILE=FileName,STATUS='Unknown')
      DO i=1,projector % numberofrows
        IF( ASSOCIATED( IntInvPerm ) ) THEN
          ii = intinvperm(i)
          IF( ii == 0 ) CYCLE
        ELSE
          ii = i
        END IF
        rowsum = 0.0_dp
        dia = 0.0_dp

        DO j=projector % rows(i), projector % rows(i+1)-1          
          jj = projector % cols(j)
          val = projector % values(j)
          IF( ii == jj ) THEN
            dia = val
          END IF
          rowsum = rowsum + val
        END DO

        IF( GlobalInds ) THEN
          ii = GlobalDofs(ii)
          WRITE(1,*) ii, i, &
              projector % rows(i+1)-projector % rows(i), ParEnv % MyPe, dia, rowsum
        ELSE
          WRITE(1,*) ii, i, &
              projector % rows(i+1)-projector % rows(i),dia, rowsum
        END IF

      END DO
      CLOSE(1)     
    END IF

    IF( ASSOCIATED(projector % rhs) ) THEN
      IF(ParEnv % PEs == 1 ) THEN
        FileName = TRIM(Prefix)//'_rhs.dat'
      ELSE
        FileName = TRIM(Prefix)//'_rhs_part'//&
            TRIM(I2S(ParEnv % MyPe))//'.dat'
      END IF
      
      OPEN(1,FILE=FileName,STATUS='Unknown')
      DO i=1,projector % numberofrows
        IF( ASSOCIATED( IntInvPerm ) ) THEN
          ii = intinvperm(i)
          IF( ii == 0 ) CYCLE
        ELSE
          ii = i
        END IF

        IF( GlobalInds ) THEN
          ii = GlobalDofs(ii)
          WRITE(1,*) ii, i, ParEnv % MyPe, projector % rhs(i)
        ELSE
          WRITE(1,*) ii, i, projector % rhs(i)
        END IF
      END DO
      CLOSE(1)     
    END IF

  END SUBROUTINE SaveProjector



  ! Set projector abs(rowsum) to unity
  !--------------------------------------------------------
  SUBROUTINE SetProjectorRowsum( Projector )
    TYPE(Matrix_t), POINTER :: Projector

    INTEGER :: i,j
    REAL(KIND=dp) :: rowsum

    DO i=1,projector % numberofrows
      rowsum = 0.0_dp
      DO j=projector % rows(i), projector % rows(i+1)-1
        rowsum = rowsum + ABS( projector % values(j) )
      END DO
      DO j=projector % rows(i), projector % rows(i+1)-1
        projector % values(j) = projector % values(j) / rowsum
      END DO
    END DO

  END SUBROUTINE SetProjectorRowsum

  
  ! This creates a projector that integrates over the BCs on the boundary such that
  ! an integral constraint may be applied on it. For example, we could set the
  ! incoming flow without actually setting the profile.
  !--------------------------------------------------------------------------------------
  FUNCTION IntegralProjector(Model, Mesh, BCInd ) RESULT ( Projector )

    TYPE(Model_t) :: Model  
    TYPE(Mesh_t), TARGET :: Mesh
    INTEGER :: BCInd    
    TYPE(Matrix_t), POINTER :: Projector
        
    REAL(KIND=dp) :: area
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found
    INTEGER :: n
    CHARACTER(*), PARAMETER :: Caller="IntegralProjector"

    
    BC => Model % BCs(BCInd) % Values
    NULLIFY(Projector)
    
    IF( .NOT. ListGetLogical( BC,'Integral BC', Found ) ) RETURN

    CALL Info(Caller,'Creating integral constraint matrix for boundary: '//TRIM(I2S(BCind)),Level=6)
    
    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_INTEGRAL
    
    CALL CreateIntegralProjector()
    
    CALL List_toCRSMatrix(Projector)
    area = SUM( Projector % Values )
    n = SIZE( Projector % Values ) 
    
    WRITE( Message,'(A,ES12.4)') 'Total area of boundary integral:',area  
    CALL Info(Caller, Message, Level=6 )

    CALL SetInvPermIndex()
    
    IF( InfoActive(20) ) THEN
       WRITE(Message,'(A,ES12.3)') 'Sum of constraint matrix entries: ',SUM(Projector%Values)
       CALL Info(Caller,Message)
       CALL Info(Caller,'Constraint matrix cols min: '//TRIM(I2S(MINVAL(Projector%Cols))))
       CALL Info(Caller,'Constraint matrix cols max: '//TRIM(I2S(MAXVAL(Projector%Cols))))
       CALL Info(Caller,'Constraint matrix rows min: '//TRIM(I2S(MINVAL(Projector%Rows))))
       CALL Info(Caller,'Constraint matrix rows max: '//TRIM(I2S(MINVAL(Projector%Rows))))
     END IF
            
  CONTAINS
    
    SUBROUTINE CreateIntegralProjector()
    
      INTEGER :: i,j,n,t,p
      REAL(KIND=dp) :: u,v,w,weight,x,detJ,val
      REAL(KIND=dp), ALLOCATABLE :: Basis(:)
      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Element
      INTEGER, POINTER :: Indexes(:)  
      TYPE(GaussIntegrationPoints_t) :: IP
      LOGICAL :: AxisSym, Stat, Visited = .FALSE.

      SAVE Visited, Nodes, Basis

      IF(.NOT. Visited ) THEN
        n = Mesh % MaxElementNodes
        ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )
        Visited = .TRUE.
      END IF

      AxisSym = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
          CurrentCoordinateSystem() == CylindricSymmetric ) 

      DO t = 1, Mesh % NumberOfBoundaryElements

        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )

        IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BCInd) % Tag ) CYCLE

        n = Element % TYPE % NumberOfNodes        
        Indexes => Element % NodeIndexes      
        IP = GaussPoints( Element )

        Nodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
        Nodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))

        DO j=1,IP % n
          u = IP % u(j)
          v = IP % v(j)
          w = IP % w(j)

          Stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis)

          weight = detJ * IP % s(j)
          IF( AxisSym ) THEN
            x = SUM( Basis(1:n) * Nodes % x(1:n) )
            weight = weight * x
          END IF
          
          DO p=1,n
            val = weight * Basis(p)
            CALL List_AddToMatrixElement(Projector % ListMatrix, 1, Indexes(p), val ) 
          END DO
          
        END DO
      END DO

    END SUBROUTINE CreateIntegralProjector    


    ! Let us associate the inverso permutation to some degree of freedom that is unique and not
    ! set by some other BCs. This unique index is needed in the future. 
    !------------------------------------------------------------------------------------------
    SUBROUTINE SetInvPermIndex()
    
      INTEGER :: i,j,t,n,maxind
      TYPE(Element_t), POINTER :: Element
      INTEGER, POINTER :: Indexes(:)  
      LOGICAL, ALLOCATABLE :: SomeOtherBC(:)
      
      IF(.NOT. ASSOCIATED( Projector % InvPerm ) ) THEN
        ALLOCATE( Projector % InvPerm(1) ) 
        Projector % InvPerm = 0
      END IF

      n = Mesh % NumberOfNodes
      ALLOCATE( SomeOtherBC(n) )
      SomeOtherBC = .FALSE.
      maxind = 0
      
      DO t = 1, Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(BCInd) % Tag ) CYCLE
        Indexes => Element % NodeIndexes      
        SomeOtherBC(Indexes) = .TRUE.
      END DO

      DO t = 1, Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(BCInd) % Tag ) THEN
          Indexes => Element % NodeIndexes      
          n = Element % TYPE % NumberOfNodes        
          DO i=1,n
            j = Indexes(i)
            IF( SomeOtherBC(j) ) CYCLE
            maxind = MAX(maxind,j)
          END DO
        END IF
      END DO
      
      IF( maxind == 0 ) THEN
        CALL Fatal(Caller,'Could not determine maximum unset index!')
      ELSE
        CALL Info(Caller,'Setting the representative node to: '//TRIM(I2S(maxind)),Level=8)
        Projector % InvPerm(1) = maxind
      END IF        
    END SUBROUTINE SetInvPermIndex
    
  END FUNCTION IntegralProjector


  
!------------------------------------------------------------------------------
!> Create a projector between Master and Target boundaries.
!> The projector may be a nodal projector x=Px or a weigted 
!> Galerking projector such that Qx=Px. In the first case the projector 
!> will be P and in the second case [Q-P]. 
!------------------------------------------------------------------------------
  FUNCTION PeriodicProjector( Model, Mesh, This, Trgt, cdim, &
      Galerkin ) RESULT(Projector)
!------------------------------------------------------------------------------   
    TYPE(Model_t) :: Model
    INTEGER :: This, Trgt
    INTEGER, OPTIONAL :: cdim
    TYPE(Mesh_t), TARGET :: Mesh
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL, OPTIONAL :: Galerkin
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,dim
    LOGICAL :: GotIt, UseQuadrantTree, Success, WeakProjector, &
        Rotational, AntiRotational, Sliding, AntiSliding, Repeating, AntiRepeating, &
        Discontinuous, NodalJump, Radial, AntiRadial, DoNodes, DoEdges, Axial, AntiAxial, &
        Flat, Plane, AntiPlane, LevelProj, FullCircle, Cylindrical, &
        ParallelNumbering, TimestepNumbering, EnforceOverlay, NormalProj
    LOGICAL, ALLOCATABLE :: MirrorNode(:)
    TYPE(Mesh_t), POINTER ::  BMesh1, BMesh2, PMesh
    TYPE(Nodes_t), POINTER :: MeshNodes, GaussNodes
    REAL(KIND=dp) :: NodeScale, EdgeScale, Radius, Coeff, val 
    TYPE(ValueList_t), POINTER :: BC
    CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix
    TYPE(Variable_t), POINTER :: v
    CHARACTER(*), PARAMETER :: Caller="PeriodicProjector"
    
    INTERFACE
      FUNCTION WeightedProjector(BMesh2, BMesh1, InvPerm2, InvPerm1, &
          UseQuadrantTree, Repeating, AntiRepeating, PeriodicScale, &
          NodalJump ) &
         RESULT ( Projector )
        USE Types
        TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
        REAL(KIND=dp) :: PeriodicScale
        INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
        LOGICAL :: UseQuadrantTree, Repeating, AntiRepeating
        TYPE(Matrix_t), POINTER :: Projector
        LOGICAL :: NodalJump
      END FUNCTION WeightedProjector
    END INTERFACE
!------------------------------------------------------------------------------
    Projector => NULL()
    IF ( This <= 0  ) RETURN    
    CALL Info(Caller,'Starting projector creation',Level=12)

    DIM = CoordinateSystemDimension()

    CALL ResetTimer(Caller)
    
    Projector => NULL()
    BC => Model % BCs(This) % Values
    PMesh => Mesh

    
    ! Whether to choose nodal or Galerkin projector is determined by an optional
    ! flag. The default is the nodal projector.
    !--------------------------------------------------------------------------
    IF( PRESENT( Galerkin) ) THEN
      WeakProjector = Galerkin
    ELSE
      WeakProjector = ListGetLogical( BC, 'Galerkin Projector', GotIt )
    END IF


    ! If the boundary is discontinuous then we have the luxury of creating the projector
    ! very cheaply using the permutation vector. This does not need the target as the 
    ! boundary is self-contained.
    !------------------------------------------------------------------------------------
    IF( ListGetLogical( BC, 'Discontinuous Boundary', GotIt ) .AND. Mesh % DisContMesh )THEN
      IF( WeakProjector ) THEN
        Projector => WeightedProjectorDiscont( PMesh, This )
      ELSE
        Projector => NodalProjectorDiscont( PMesh, This )
      END IF
      
      IF ( .NOT. ASSOCIATED( Projector ) ) RETURN
      GOTO 100
    END IF
    
    IF ( Trgt <= 0 ) RETURN    

    ! Create the mesh projector, and if needed, also eliminate the ghost nodes
    ! There are two choices of projector: a nodal projector P in x=Px, and a 
    ! Galerkin projector [Q-P] in Qx=Px. 
    ! The projector is assumed to be either a rotational projector with no translation
    ! and rotation, or then generic one with possible coordinate mapping.
    !---------------------------------------------------------------------------------
    CALL Info(Caller,'-----------------------------------------------------',Level=8)
    WRITE( Message,'(A,I0,A,I0)') 'Creating projector between BCs ',This,' and ',Trgt
    CALL Info(Caller,Message,Level=8)

    ! Create temporal mesh structures that are utilized when making the 
    ! projector between "This" and "Trgt" boundary.
    !--------------------------------------------------------------------------
    BMesh1 => AllocateMesh()
    BMesh2 => AllocateMesh()

    CALL CreateInterfaceMeshes( Model, Mesh, This, Trgt, Bmesh1, BMesh2, &
        Success ) 

    IF(.NOT. Success) THEN
      CALL Info(Caller,'Releasing interface meshes!',Level=20)
      CALL ReleaseMesh(BMesh1)
      CALL ReleaseMesh(BMesh2)
      RETURN
    END IF

    ! If requested map the interface coordinate from (x,y,z) to any permutation of these. 
    CALL MapInterfaceCoordinate( BMesh1, BMesh2, Model % BCs(This) % Values )

    NormalProj = ListGetLogical( BC,'Normal Projector',GotIt )
    
    ! Check whether to use (anti)rotational projector.
    ! We don't really know on which side the projector was called so 
    ! let's check both sides.
    !--------------------------------------------------------------------------
    Rotational = ListGetLogical( BC,'Rotational Projector',GotIt )
    AntiRotational = ListGetLogical( BC,'Anti Rotational Projector',GotIt )
    IF( AntiRotational ) Rotational = .TRUE.

    Cylindrical =  ListGetLogical( BC,'Cylindrical Projector',GotIt )

    Radial = ListGetLogical( BC,'Radial Projector',GotIt )
    AntiRadial = ListGetLogical( BC,'Anti Radial Projector',GotIt )
    IF( AntiRadial ) Radial = .TRUE.

    Axial = ListGetLogical( BC,'Axial Projector',GotIt )
    AntiAxial = ListGetLogical( BC,'Anti Axial Projector',GotIt )
    IF( AntiAxial ) Axial = .TRUE.

    Sliding = ListGetLogical( BC,'Sliding Projector',GotIt )
    AntiSliding = ListGetLogical( BC,'Anti Sliding Projector',GotIt )
    IF( AntiSliding ) Sliding = .TRUE. 

    Flat = ListGetLogical( BC,'Flat Projector',GotIt )
    Plane = ListGetLogical( BC, 'Plane Projector',GotIt )
    AntiPlane = ListGetLogical( BC,'Anti Plane Projector',GotIt )    
    IF( AntiPlane ) Plane = .TRUE.
    
    IF( Radial ) CALL Info(Caller,'Enforcing > Radial Projector <',Level=12)
    IF( Axial ) CALL Info(Caller,'Enforcing > Axial Projector <',Level=12)
    IF( Sliding ) CALL Info(Caller,'Enforcing > Sliding Projector <',Level=12)
    IF( Cylindrical ) CALL Info(Caller,'Enforcing > Cylindrical Projector <',Level=12)
    IF( Rotational ) CALL Info(Caller,'Enforcing > Rotational Projector <',Level=12)
    IF( Flat ) CALL Info(Caller,'Enforcing > Flat Projector <',Level=12)
    IF( Plane ) CALL Info(Caller,'Enforcing > Plane Projector <',Level=12)

    NodeScale = ListGetConstReal( BC, 'Mortar BC Scaling',GotIt)
    IF(.NOT.Gotit ) THEN
      IF( AntiRadial .OR. AntiPlane ) THEN
        NodeScale = -1._dp
      ELSE
        NodeScale = 1.0_dp
      END IF
    END IF
    EdgeScale = NodeScale

    NodalJump = ListCheckPrefix( BC,'Mortar BC Coefficient')
    IF(.NOT. NodalJump ) THEN
      NodalJump = ListCheckPrefix( BC,'Mortar BC Resistivity')
    END IF

    ! There are tailored projectors for simplified interfaces
    !-------------------------------------------------------------

    ! Stride projector is obsolete and has been eliminated.
    IF( ListGetLogical( BC,'Stride Projector',GotIt) ) THEN
      CALL ListAddLogical( BC,'Level Projector',.TRUE.)
      CALL ListAddLogical( BC,'Level Projector Strong',.TRUE.)
      CALL Warn(Caller,'Enforcing > Level Projector < instead of old > Stride Projector <')
    END IF

    LevelProj = ListGetLogical( BC,'Level Projector',GotIt) 
    IF( Rotational .OR. Cylindrical .OR. Radial .OR. Flat .OR. Plane .OR. Axial ) THEN
      IF(.NOT. GotIt ) THEN
        CALL Info(Caller,'Enforcing > Level Projector = True < with dimensional reduction',&
            Level = 7 )
        LevelProj = .TRUE. 
      ELSE IF(.NOT. LevelProj ) THEN
        ! If we have dimensionally reduced projector but don't use LevelProjector 
        ! to integrate over it, then ensure that the 3rd coordinate is set to zero.
        BMesh1 % Nodes % z = 0.0_dp
        BMesh2 % Nodes % z = 0.0_dp
      END IF
    END IF


    IF( LevelProj ) THEN
      IF( ListGetLogical( Model % Solver % Values,'Projector Skip Nodes',GotIt ) ) THEN
        DoNodes = .FALSE.
      ELSE
        IF( ListGetLogical( BC,'Projector Skip Nodes',GotIt) ) THEN
          DoNodes = .FALSE.
        ELSE
          DoNodes = ( Mesh % NumberOfNodes > 0 ) 
        END IF
      END IF

      IF( ListGetLogical( Model % Solver % Values,'Projector Skip Edges',GotIt ) ) THEN
        DoEdges = .FALSE.
      ELSE
        IF( ListGetLogical( BC,'Projector Skip Edges',GotIt) ) THEN
          DoEdges = .FALSE.
        ELSE
          ! We are conservative here since there may be edges in 2D which 
          ! still cannot be used for creating the projector
          DoEdges = ( Mesh % NumberOfEdges > 0 .AND. &
              Mesh % MeshDim == 3 .AND. Dim == 3 )

          ! Ensure that there is no p-elements that made us think that we have edges
          ! Here we assume that if there is any p-element then also the 1st element is such
          IF( DoEdges ) THEN
            IF(isPelement(Mesh % Elements(1))) THEN
              DoEdges = .FALSE.
              CALL Info(Caller,'Edge projector will not be created for p-element mesh',Level=10)
            END IF
          END IF
        END IF
      END IF
    END IF


    ! If the interface is rotational move to (phi,z) plane and alter the phi coordinate
    ! so that the meshes coincide.
    ! Otherwise make the two meshes to coincide using rotation, translation &
    ! scaling.
    !---------------------------------------------------------------------------------
    Radius = 1.0_dp
    FullCircle = .FALSE.
    EnforceOverlay = ListGetLogical( BC, 'Mortar BC enforce overlay', GotIt )

    IF( .NOT. Rotational ) THEN
      IF( ListCheckPresent( BC,'Mesh Rotate 3' ) ) THEN
        CALL Fatal(Caller,'Only "Rotational Projector" has the "Mesh Rotate 3" trick implemented')
      END IF
    END IF

    IF( Rotational .OR. Cylindrical ) THEN
      CALL RotationalInterfaceMeshes( BMesh1, BMesh2, BC, Cylindrical, &
          Radius, FullCircle )
    ELSE IF( Radial ) THEN
      CALL RadialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Flat ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Axial ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )      
      CALL AxialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Plane ) THEN
      CALL PlaneInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( .NOT. ( Sliding .OR. NormalProj ) ) THEN
      IF( .NOT. GotIt ) EnforceOverlay = .TRUE.
    END IF

    IF( EnforceOverlay ) THEN
      CALL OverlayIntefaceMeshes( BMesh1, BMesh2, BC )
    END IF

    Repeating = ( Rotational .OR. Sliding .OR. Axial ) .AND. .NOT. FullCircle 
    AntiRepeating = .FALSE.
    IF( Repeating ) THEN
      AntiRepeating = ListGetLogical( BC,'Antisymmetric BC',GotIt ) 
      IF( .NOT. GotIt ) THEN
        AntiRepeating = ( AntiRotational .OR. AntiSliding .OR. AntiAxial ) .AND. .NOT. FullCircle 
      END IF
    END IF
      
    IF( LevelProj ) THEN 
      Projector => LevelProjector( BMesh1, BMesh2, Repeating, AntiRepeating, &
          FullCircle, Radius, DoNodes, DoEdges, &          
          NodeScale, EdgeScale, BC )
    ELSE IF( NormalProj ) THEN
      IF( AntiRepeating ) THEN
        CALL Fatal(Caller,'An antiperiodic projector cannot be dealt with the normal projector!')
      END IF
      Projector => NormalProjector( BMesh2, BMesh1, BC )
    ELSE
      IF( FullCircle ) THEN
        CALL Fatal(Caller,'A full circle cannot be dealt with the generic projector!')
      END IF

      UseQuadrantTree = ListGetLogical(Model % Simulation,'Use Quadrant Tree',GotIt)
      IF( .NOT. GotIt ) UseQuadrantTree = .TRUE.
      
      IF( WeakProjector ) THEN
        CALL Info(Caller,'Weighted projector is very suboptimal, do not use it!',Level=3)
        CALL Warn(Caller,'Use "Normal Projector" or "Level Projector" instead!')
        Projector => WeightedProjector( BMesh2, BMesh1, BMesh2 % InvPerm, BMesh1 % InvPerm, &
            UseQuadrantTree, Repeating, AntiRepeating, NodeScale, NodalJump )
      ELSE
        Projector => NodalProjector( BMesh2, BMesh1, &
            UseQuadrantTree, Repeating, AntiRepeating )
      END IF
    END IF

    IF( InfoActive(15) ) THEN
      val = SUM( Projector % Values )
      WRITE(Message,'(A,ES12.3)') 'Sum of projector entries:',val
      CALL Info(Caller,Message)
      
      val = MINVAL( Projector % Values )
      WRITE(Message,'(A,ES12.3)') 'Minimum of projector entries:',val
      CALL Info(Caller,Message)
      
      val = MAXVAL( Projector % Values )
      WRITE(Message,'(A,ES12.3)') 'Maximum of projector entries:',val
      CALL Info(Caller,Message)
      
      CALL Info(Caller,'Number of rows in projector: '&
          //TRIM(I2S(Projector % NumberOfRows)))
      CALL Info(Caller,'Number of entries in projector: '&
          //TRIM(I2S(SIZE(Projector % Values))))
    END IF


    
    ! Deallocate mesh structures:
    !---------------------------------------------------------------
    CALL Info(Caller,'Releasing interface meshes!',Level=20)
    BMesh1 % Projector => NULL()
    BMesh1 % Parent => NULL()
    !DEALLOCATE( BMesh1 % InvPerm ) 
    CALL ReleaseMesh(BMesh1)

    BMesh2 % Projector => NULL()
    BMesh2 % Parent => NULL()
    !DEALLOCATE( BMesh2 % InvPerm ) 
    CALL ReleaseMesh(BMesh2)

100 Projector % ProjectorBC = This

    IF( ListGetLogical( BC,'Projector Set Rowsum',GotIt ) ) THEN
      CALL SetProjectorRowsum( Projector )
    END IF

    Coeff = ListGetConstReal( BC,'Projector Multiplier',GotIt) 
    IF(.NOT. GotIt) Coeff = ListGetConstReal( Model % Simulation,&
        'Projector Multiplier',GotIt) 
    IF( GotIt ) Projector % Values = Coeff * Projector % Values

    IF( ListGetLogical( BC,'Save Projector',GotIt ) ) THEN
      ParallelNumbering = ListGetLogical( BC,'Save Projector Global Numbering',GotIt )

      FilePrefix = 'p'//TRIM(I2S(This))
      
      TimestepNumbering = ListGetLogical( BC,'Save Projector Timestep Numbering',GotIt )
      IF( TimestepNumberIng ) THEN
        i = 0
        v => VariableGet( Mesh % Variables, 'timestep' )
        IF( ASSOCIATED( v ) ) i = NINT( v % Values(1) )
        WRITE( FilePrefix,'(A,I4.4)') TRIM(FilePrefix)//'_',i
      END IF
        
      CALL SaveProjector( Projector, .TRUE.,TRIM(FilePrefix), &
          Parallel = ParallelNumbering) 
      
      ! Dual projector if it exists
      IF( ASSOCIATED( Projector % Ematrix ) ) THEN
        CALL SaveProjector( Projector % Ematrix, .TRUE.,'dual_'//TRIM(FilePrefix),&
            Projector % InvPerm, Parallel = ParallelNumbering) 
      END IF

      ! Biorthogonal projector if it exists
      IF( ASSOCIATED( Projector % Child ) ) THEN
        CALL SaveProjector( Projector % Child, .TRUE.,'biortho_'//TRIM(FilePrefix), &
            Projector % InvPerm, Parallel = ParallelNumbering ) 
      END IF

      IF( ListGetLogical( BC,'Save Projector And Stop',GotIt ) ) STOP EXIT_OK
    END IF    

    CALL CheckTimer(Caller,Delete=.TRUE.)
    CALL Info(Caller,'Projector created, now exiting...',Level=8)

!------------------------------------------------------------------------------
  END FUNCTION PeriodicProjector
!------------------------------------------------------------------------------


  

!------------------------------------------------------------------------------
!> Create a permutation between two meshes such that we can solve a smaller system.
!------------------------------------------------------------------------------
  SUBROUTINE PeriodicPermutation( Model, Mesh, This, Trgt, PerPerm, PerFlip, DoFaces ) 
!------------------------------------------------------------------------------   
    TYPE(Model_t) :: Model
    INTEGER :: This, Trgt
    TYPE(Mesh_t), TARGET :: Mesh
    INTEGER, POINTER :: PerPerm(:)
    LOGICAL, POINTER :: PerFlip(:)
    LOGICAL :: DoFaces
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,dim
    LOGICAL :: GotIt, Success, Rotational, AntiRotational, Sliding, AntiSliding, Repeating, &
        Radial, AntiRadial, DoNodes, DoEdges, Axial, AntiAxial, &
        Flat, Plane, AntiPlane, Cylindrical, ParallelNumbering, EnforceOverlay, &
        FullCircle, AntiPeriodic
    REAL(KIND=dp) :: Radius
    TYPE(Mesh_t), POINTER ::  BMesh1, BMesh2, PMesh
    TYPE(ValueList_t), POINTER :: BC
    
!------------------------------------------------------------------------------
    IF ( This <= 0  .OR. Trgt <= 0 ) RETURN    
    CALL Info('PeriodicPermutation','Starting periodic permutation creation',Level=12)

    CALL ResetTimer('PeriodicPermutation')
    
    DIM = CoordinateSystemDimension()
    BC => Model % BCs(This) % Values
    PMesh => Mesh
    
    CALL Info('PeriodicPermutation','-----------------------------------------------------',Level=8)
    WRITE( Message,'(A,I0,A,I0)') 'Creating mapping between BCs ',This,' and ',Trgt
    CALL Info('PeriodicPermutation',Message,Level=8)

    BMesh1 => AllocateMesh()
    BMesh2 => AllocateMesh()
    
    CALL CreateInterfaceMeshes( Model, Mesh, This, Trgt, Bmesh1, BMesh2, Success ) 
    
    IF(.NOT. Success) THEN
      CALL Info('PeriodicPermutation','Releasing interface meshes!',Level=20)
      CALL ReleaseMesh(BMesh1)
      CALL ReleaseMesh(BMesh2)
      RETURN
    END IF

    ! If requested map the interface coordinate from (x,y,z) to any permutation of these. 
    CALL MapInterfaceCoordinate( BMesh1, BMesh2, Model % BCs(This) % Values )
    
    ! Lets check what kind of symmetry we have.
    Rotational = ListGetLogical( BC,'Rotational Projector',GotIt )
    AntiRotational = ListGetLogical( BC,'Anti Rotational Projector',GotIt )

    Cylindrical =  ListGetLogical( BC,'Cylindrical Projector',GotIt )

    Radial = ListGetLogical( BC,'Radial Projector',GotIt )
    AntiRadial = ListGetLogical( BC,'Anti Radial Projector',GotIt )
    IF( AntiRadial ) Radial = .TRUE.
    
    Axial = ListGetLogical( BC,'Axial Projector',GotIt )
    AntiAxial = ListGetLogical( BC,'Anti Axial Projector',GotIt )
    IF( AntiAxial ) Axial = .TRUE.
    
    Sliding = ListGetLogical( BC, 'Sliding Projector',GotIt )
    AntiSliding = ListGetLogical( BC, 'Anti Sliding Projector',GotIt )
    IF( AntiSliding ) Sliding = .TRUE.
    
    Flat = ListGetLogical( BC, 'Flat Projector',GotIt )
    Plane = ListGetLogical( BC, 'Plane Projector',GotIt )
    AntiPlane = ListGetLogical( BC,'Anti Plane Projector',GotIt )    
    IF( AntiPlane ) Plane = .TRUE.

    AntiPeriodic = ListGetLogical( BC,'Antisymmetric BC',GotIt )
    IF( .NOT. GotIt ) THEN   
      AntiPeriodic = ( AntiRotational .OR. AntiRadial .OR. AntiAxial .OR. AntiPlane ) 
    END IF
      
    IF( AntiPeriodic ) CALL Info('PeriodicPermutation','Assuming antiperiodic conforming projector',Level=8)
    
    IF( Radial ) CALL Info('PeriodicPermutation','Enforcing > Radial Projector <',Level=12)
    IF( Axial ) CALL Info('PeriodicPermutation','Enforcing > Axial Projector <',Level=12)
    IF( Sliding ) CALL Info('PeriodicPermutation','Enforcing > Sliding Projector <',Level=12)
    IF( Cylindrical ) CALL Info('PeriodicPermutation','Enforcing > Cylindrical Projector <',Level=12)
    IF( Rotational ) CALL Info('PeriodicPermutation','Enforcing > Rotational Projector <',Level=12)
    IF( Flat ) CALL Info('PeriodicPermutation','Enforcing > Flat Projector <',Level=12)
    IF( Plane ) CALL Info('PeriodicPermutation','Enforcing > Plane Projector <',Level=12)

    DoNodes = .TRUE.
    !IF( ListGetLogical( Model % Solver % Values,'Projector Skip Nodes',GotIt ) ) DoNodes = .FALSE.    
    IF( ListGetLogical( BC,'Projector Skip Nodes',GotIt) ) DoNodes = .FALSE.

    ! We are conservative here since there may be edges in 2D which 
    ! still cannot be used for creating the projector
    DoEdges = ( Mesh % NumberOfEdges > 0 .AND. Mesh % MeshDim == 3 .AND. Dim == 3 )
    
    ! Ensure that there is no p-elements that made us think that we have edges
    ! Here we assume that if there is any p-element then also the 1st element is such
    IF( DoEdges ) THEN
      IF(isPelement(Mesh % Elements(1))) THEN
        DoEdges = .FALSE.
        CALL Info('PeriodicPermutation','Edge projector will not be created for p-element mesh',Level=10)
      END IF
    END IF
        
    !IF( ListGetLogical( Model % Solver % Values,'Projector Skip Edges',GotIt ) ) DoEdges = .FALSE.
    IF( ListGetLogical( BC,'Projector Skip Edges',GotIt) ) DoEdges = .FALSE.
      
    ! Make the two meshes to coincide using rotation, translation scaling.
    !---------------------------------------------------------------------------------
    Radius = 1.0_dp
    EnforceOverlay = ListGetLogical( BC, 'Mortar BC enforce overlay', GotIt )

    IF( Rotational .OR. Cylindrical ) THEN
      CALL RotationalInterfaceMeshes( BMesh1, BMesh2, BC, Cylindrical, &
          Radius, FullCircle )
      IF( FullCircle ) CALL Fatal('PeriodicPermutation','Cannot deal full circle with permutation')
    ELSE IF( Radial ) THEN
      CALL RadialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Flat ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Axial ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )      
      CALL AxialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Plane ) THEN
      CALL PlaneInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( .NOT. Sliding ) THEN
      IF( .NOT. GotIt ) EnforceOverlay = .TRUE.
    END IF
    
    IF( EnforceOverlay ) THEN
      CALL OverlayIntefaceMeshes( BMesh1, BMesh2, BC )
    END IF
    
    IF( DoNodes ) CALL ConformingNodePerm(PMesh, BMesh1, BMesh2, PerPerm, PerFlip, AntiPeriodic )
    IF( DoEdges ) CALL ConformingEdgePerm(PMesh, BMesh1, BMesh2, PerPerm, PerFlip, AntiPeriodic )
    IF( DoEdges .AND. DoFaces ) &
        CALL ConformingFacePerm(PMesh, BMesh1, BMesh2, PerPerm, PerFlip, AntiPeriodic )
    
    ! Deallocate mesh structures:
    !---------------------------------------------------------------
    CALL Info('PeriodicPermutation','Releasing interface meshes!',Level=20)
    BMesh1 % Projector => NULL()
    BMesh1 % Parent => NULL()
    !DEALLOCATE( BMesh1 % InvPerm ) 
    CALL ReleaseMesh(BMesh1)

    BMesh2 % Projector => NULL()
    BMesh2 % Parent => NULL()
    !DEALLOCATE( BMesh2 % InvPerm ) 
    CALL ReleaseMesh(BMesh2)

    CALL CheckTimer('PeriodicPermutation',Delete=.TRUE.)
           
    CALL Info('PeriodicPermutation','Periodic permutation created, now exiting...',Level=8)
   
    
!------------------------------------------------------------------------------
  END SUBROUTINE PeriodicPermutation
!------------------------------------------------------------------------------


  
  !> If periodic BCs given, compute boundary mesh projector.
  !> If conforming BCs given, create permutation for elimination.
  !------------------------------------------------------
  SUBROUTINE GeneratePeriodicProjectors( Model, Mesh ) 
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: i,j,k,n,nocyclic,noconf,noflip,mini,maxi
    LOGICAL :: Found, NeedFaces
    INTEGER, POINTER :: PerPerm(:)
    LOGICAL, POINTER :: PerFlip(:)

    TYPE(Mesh_t), POINTER :: sMesh
    TYPE(Solver_t), POINTER :: sSolver
 
! set these to satisfy possible call to EdgeElementInfo() - and restore later. Maybe not very
! beautiful, but seems to work for now.
! x
    sSolver => Model % Solver
    Model % Solver => Model % Solvers(1)
    sMesh => Model % Solver % Mesh
    Model % Solver % Mesh => Mesh
! x

    
    DO i = 1,Model % NumberOfBCs
      k = ListGetInteger( Model % BCs(i) % Values, 'Periodic BC', Found )
      IF( Found ) THEN
        Model % BCs(i) % PMatrix => PeriodicProjector( Model, Mesh, i, k )
      END IF
    END DO

    
    IF( ListCheckPresentAnyBC( Model,'Conforming BC' ) ) THEN
      NeedFaces = ListGetLogicalAnySolver( Model,'Use Piola Transform')
      
      n = Mesh % NumberOfNodes + Mesh % NumberOfEdges
      IF(NeedFaces) n = n + 2 * Mesh % NumberOfFaces + 3 * Mesh % NumberOfBulkElements

      IF(.NOT. ASSOCIATED( Mesh % PeriodicPerm ) ) THEN
        CALL Info('GeneratePeriodicProjectors','Allocating for conforming data!')
        ALLOCATE( Mesh % PeriodicPerm(n) )
        ALLOCATE( Mesh % PeriodicFlip(n) )
      END IF      
      PerPerm => Mesh % PeriodicPerm      
      PerPerm = 0
      PerFlip => Mesh % PeriodicFlip
      PerFlip = .FALSE.

      DO i = 1,Model % NumberOfBCs
        k = ListGetInteger( Model % BCs(i) % Values, 'Conforming BC', Found )
        IF( Found ) THEN
          CALL PeriodicPermutation( Model, Mesh, i, k, PerPerm, PerFlip, NeedFaces )
        END IF
      END DO
      nocyclic = 0
      noconf = 0
      mini = HUGE(mini)
      maxi = 0
      
      DO i = 1,n
        j = PerPerm(i)
        IF( j > 0 ) THEN
          mini = MIN( mini, i )
          maxi = MAX( maxi, i )
          noconf = noconf + 1
          IF( PerPerm(j) > 0 ) THEN
            PerPerm(i) = PerPerm(j)
            IF( PerFlip(i) ) THEN
              PerFlip(i) = .NOT. PerFlip(j)
            ELSE
              PerFlip(i) = PerFlip(j)
            END IF
            nocyclic = nocyclic + 1
          END IF
        END IF
      END DO
      noflip = COUNT( PerFlip )
            
      CALL Info('GeneratePeriodicProjectors','Number of conforming maps: '//TRIM(I2S(noconf)),Level=8)
      IF(nocyclic>0) CALL Info('GeneratePeriodicProjectors','Number of cyclic maps: '//TRIM(I2S(nocyclic)),Level=8)
      IF(noflip>0) CALL Info('GeneratePeriodicProjectors','Number of periodic flips: '//TRIM(I2S(noflip)),Level=8)
    END IF

    Model % Solver % Mesh => sMesh
    Model % Solver => sSolver

    
  END SUBROUTINE GeneratePeriodicProjectors


!------------------------------------------------------------------------------
!> Create node distribution for a unit segment x \in [0,1] with n elements 
!> i.e. n+1 nodes. There are different options for the type of distribution.
!> 1) Even distribution 
!> 2) Geometric distribution
!> 3) Arbitrary distribution determined by a functional dependence
!> Note that the 3rd algorithm involves iterative solution of the nodal
!> positions and is therefore not bullet-proof.
!------------------------------------------------------------------------------
  SUBROUTINE UnitSegmentDivision( w, n, ExtList )
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    INTEGER :: n
    TYPE(ValueList_t), POINTER, OPTIONAL :: ExtList
    !---------------------------------------------------------------
    INTEGER :: i,J,iter,maxiter
    REAL(KIND=dp) :: q,r,h1,hn,minhn,err_eps,err,xn
    REAL(KIND=dp), ALLOCATABLE :: wold(:),h(:)
    LOGICAL :: Found, GotRatio, FunExtruded, Fun1D
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: ParList
    
    IF( PRESENT( ExtList ) ) THEN
      ParList => ExtList
    ELSE
      ParList => CurrentModel % Simulation
    END IF

    FunExtruded = ListCheckPresent( ParList,'Extruded Mesh Density')
    Fun1D = ListCheckPresent( ParList,'1D Mesh Density')
    
    ! Geometric division
    !---------------------------------------------------------------
    q = ListGetConstReal( ParList,'Extruded Mesh Ratio',GotRatio)
    IF(.NOT. GotRatio) q = ListGetConstReal( ParList,'1D Mesh Ratio',GotRatio)
    IF( GotRatio ) THEN
      IF( ( ABS(ABS(q)-1.0_dp) < 1.0d-6 ) .OR. (q < 0.0_dp .AND. n <= 2) ) THEN
        CALL Info('UnitSegmentDivision','Assuming linear division as mesh ratio is close to one!')
        GotRatio = .FALSE.
      END IF
    END IF
    
    IF( GotRatio ) THEN
      CALL Info('UnitSegmentDivision','Creating geometric division',Level=5)

      IF( q > 0.0_dp ) THEN      
        r = q**(1.0_dp/(n-1))
        h1 = (1-r)/(1-r**n)
        w(0) = 0.0_dp
        DO i=1,n-1
          w(i) = h1 * (1-r**i)/(1-r)
        END DO
        w(n) = 1.0_dp
      ELSE
        q = -q
        IF(MODULO(n,2) == 0) THEN
          r = q**(1.0_dp/(n/2-1))
          h1 = 0.5_dp*(1-r)/(1-r**(n/2))
        ELSE 
          r = q**(1.0_dp/((n-1)/2))
          h1 = 0.5_dp / ( (1-r**((n+1)/2))/(1-r) - 0.5_dp * r**((n-1)/2))
        END IF
        
        w(0) = 0.0_dp
        DO i=1,n
          IF( i <= n/2 ) THEN
            w(i) = h1 * (1-r**i)/(1-r)
          ELSE
            w(i) = 1.0_dp -  h1 * (1-r**(n-i))/(1-r)
          END IF
        END DO
        w(n) = 1.0_dp
      END IF
            
    ! Generic division given by a function
    !-----------------------------------------------------------------------
    ELSE IF( FunExtruded .OR. Fun1D ) THEN

      CALL Info('UnitSegmentDivision','Creating functional division',Level=5)

      ! Initial guess is an even distribution
      DO i=0,n
        w(i) = i/(1._dp * n)
      END DO

      ALLOCATE( wold(0:n),h(1:n))
      wold = w

      ! parameters that determine the accuracy of the iteration
      maxiter = 10000
      err_eps = 1.0d-6

      ! Iterate to have a density distribution
      !---------------------------------------
      DO iter=1,maxiter
        
        minhn = HUGE(minhn)
        wold = w

        ! Compute the point in the local mesh xn \in [0,1]  
        ! and get the mesh parameter for that element from
        ! external function.
        !---------------------------------------------------
        DO i=1,n
          xn = (w(i)+w(i-1))/2.0_dp
          minhn = MIN( minhn, w(i)-w(i-1) )
          IF( FunExtruded ) THEN
            h(i) = ListGetFun( ParList,'Extruded Mesh Density', xn )
          ELSE
            h(i) = ListGetFun( ParList,'1D Mesh Density', xn )
          END IF
          IF( h(i) < EPSILON( h(i) ) ) THEN
            CALL Fatal('UnitSegmentDivision','Given value for h(i) was negative!')
          END IF
        END DO

        ! Utilize symmetric Gauss-Seidel to compute the new positions, w(i).
        ! from a weigted mean of the desired elemental densities, h(i).
        ! Note that something more clever could be applied here. 
        ! This was just a first implementation...
        !-------------------------------------------------------------
        DO i=1,n-1
          w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
        END DO
        DO i=n-1,1,-1
          w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
        END DO
        
        ! If the maximum error is small compared to the minimum elementsize then exit
        !-----------------------------------------------------------------------------
        err = MAXVAL( ABS(w-wold))/minhn

        IF( err < err_eps ) THEN
          WRITE( Message, '(A,I0,A)') 'Convergence obtained in ',iter,' iterations'
          CALL Info('UnitSegmentDivision', Message, Level=9 )
          EXIT
        END IF
      END DO

      IF( iter > maxiter ) THEN
        CALL Warn('UnitSegmentDivision','No convergence obtained for the unit mesh division!')
      END IF

    ! Uniform division 
    !--------------------------------------------------------------
    ELSE
      CALL Info('UnitSegmentDivision','Creating linear division',Level=5)
      DO i=0,n     
        w(i) = i/(1._dp * n)
      END DO
    END IF
    
    CALL Info('UnitSegmentDivision','Mesh division ready',Level=9)
    DO i=0,n
      WRITE( Message, '(A,I0,A,ES12.4)') 'w(',i,') : ',w(i)
      CALL Info('UnitSegmentDivision', Message, Level=9 )
    END DO

  END SUBROUTINE UnitSegmentDivision
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Given a 2D mesh extrude it to be 3D. The 3rd coordinate will always
!> be at the interval [0,1]. Therefore the adaptation for different shapes
!> must be done with StructuredMeshMapper, or some similar utility. 
!> The top and bottom surface will be assigned Boundary Condition tags
!> with indexes one larger than the maximum used on by the 2D mesh. 
!> NOTE: This function handles NDOFs of the element structure in a way
!>       which is not consistent with "Element = n:N ...", with N>1 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrude(Mesh_in, in_levels) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh_in, Mesh_out
    INTEGER :: in_levels
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: ExtrudedMeshName
    INTEGER :: i,j,k,l,n,cnt,cnt101,ind(8),max_baseline_bid,max_bid,l_n,max_body,bcid,&
        ExtrudedCoord,dg_n,totalnumberofelements
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    INTEGER :: nnodes,gnodes,gelements,ierr
    LOGICAL :: isParallel, Found, NeedEdges, PreserveBaseline, PreserveEdges, &
        Rotational, Rotate2Pi
    REAL(KIND=dp)::w,MinCoord,MaxCoord,CurrCoord
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
!------------------------------------------------------------------------------

    CALL Info('MeshExtrude','Creating '//TRIM(I2S(in_levels+1))//' extruded element layers',Level=10)

    Mesh_out => AllocateMesh()

    isParallel = ParEnv % PEs>1

    ! Generate volume nodal points:
    ! -----------------------------
    n=Mesh_in % NumberOfNodes
    nnodes=(in_levels+2)*n
    gnodes = nnodes

    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )

    gelements = Mesh_in % NumberOfBulkElements

    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo
    
      IF(.NOT. ASSOCIATED( PI_in ) ) CALL Fatal('MeshExtrude','PI_in not associated!')
      IF(.NOT. ASSOCIATED( PI_out ) ) CALL Fatal('MeshExtrude','PI_out not associated!')
            
      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % NodeInterface(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      IF(.NOT. ASSOCIATED( PI_in % NeighbourList ) ) THEN
        CALL Fatal('MeshExtrude','Neighnours not associated!')
      END IF

      ! For unset neighbours just set the this partition to be the only owner
      DO i=1,Mesh_in % NumberOfNodes
        IF (.NOT.ASSOCIATED(PI_in % NeighbourList(i) % Neighbours)) THEN
          CALL AllocateVector(PI_in % NeighbourList(i) % Neighbours,1)
          PI_in % NeighbourList(i) % Neighbours(1) = ParEnv % Mype
        END IF
      END DO
          
      j=0
      DO i=1,Mesh_in % NumberOfNodes
        IF (PI_in % NeighbourList(i) % &
            Neighbours(1) == ParEnv % MyPE ) j=j+1
      END DO

      CALL MPI_ALLREDUCE(j,gnodes,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
      
      j=0
      DO i=1,Mesh_in % NumberOfBulkElements
        IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
      END DO
      
      CALL MPI_ALLREDUCE(j,gelements,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF

    CALL Info('MeshExtrude','Number of extruded nodes: '//TRIM(I2S(nnodes)),Level=12)
    CALL Info('MeshExtrude','Number of extruded elements: '//TRIM(I2S(gelements)),Level=12)


    ! Create the division for the 1D unit mesh
    !--------------------------------------------
    ALLOCATE( Wtable( 0: in_levels + 1 ) )
    CALL UnitSegmentDivision( Wtable, in_levels + 1 ) 

    ExtrudedCoord = ListGetInteger( CurrentModel % Simulation,'Extruded Coordinate Index', &
        Found, minv=1,maxv=3 )
    IF(.NOT. Found) ExtrudedCoord = 3 

    IF( ExtrudedCoord == 1 ) THEN
      ActiveCoord => Mesh_out % Nodes % x
    ELSE IF( ExtrudedCoord == 2 ) THEN
      ActiveCoord => Mesh_out % Nodes % y
    ELSE IF( ExtrudedCoord == 3 ) THEN
      ActiveCoord => Mesh_out % Nodes % z
    END IF


    PreserveBaseline = ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found )
    IF(.NOT. Found) PreserveBaseline = .FALSE.

    PreserveEdges = ListGetLogical( CurrentModel % Simulation,'Preserve Edges',Found )
    IF(.NOT. Found) PreserveEdges = .FALSE.

    MinCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Min Coordinate',Found )
    IF(.NOT. Found) MinCoord = 0.0_dp

    MaxCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Max Coordinate',Found )
    IF(.NOT. Found) MaxCoord = 1.0_dp

    Rotate2Pi = .FALSE.
    Rotational = ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Rotational',Found )    
    IF( Rotational ) THEN
      Rotate2Pi = ( ABS(ABS( MaxCoord-MinCoord ) - 2*PI) < 1.0d-3*PI )
      IF( Rotate2Pi ) CALL Info('MeshExtrude','Perfoming full 2Pi rotation',Level=6)
    END IF

    
    cnt=0
    DO i=0,in_levels+1

      ! If we rotate full 2Pi then we have natural closure!
      IF( Rotate2Pi ) THEN
        IF( i == in_levels+1) EXIT
      END IF
      
      w = Wtable( i ) 
      CurrCoord = w * MaxCoord + (1-w) * MinCoord      
      
      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          PI_out % NodeInterface(cnt) = PI_in % NodeInterface(j)

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(&
               SIZE(PI_in % NeighbourList(j) % Neighbours)))
          PI_out % NeighbourList(cnt) % Neighbours = &
            PI_in % NeighbourList(j) % Neighbours

          PI_out % GlobalDOFs(cnt) = PI_in % GlobalDOFs(j)+i*gnodes
        END IF

      END DO
    END DO
    Mesh_out % NumberOfNodes=cnt
    Mesh_out % Nodes % NumberOfNodes = cnt

    
    IF( Rotational ) THEN
      BLOCK
        REAL(KIND=DP) :: x,y,z,r        
        DO i=1,cnt          
          x = Mesh_out % Nodes % x(i)
          y = Mesh_out % Nodes % y(i)
          z = Mesh_out % Nodes % z(i)

          Mesh_out % Nodes % x(i) = COS(z) * x
          Mesh_out % Nodes % y(i) = SIN(z) * x
          Mesh_out % Nodes % z(i) = y
        END DO
      END BLOCK
    END IF
    
    
    ! Count 101 elements:
    ! (these require an extra layer)
    ! -------------------

    cnt101 = 0
    DO i=Mesh_in % NumberOfBulkElements+1, &
         Mesh_in % NumberOfBulkElements+Mesh_in % NumberOfBoundaryElements
       IF(Mesh_in % Elements(i) % TYPE % ElementCode == 101) cnt101 = cnt101+1
    END DO

    n=SIZE(Mesh_in % Elements)

    ! inquire total number of needed 
    IF( Rotate2Pi ) THEN
      totalnumberofelements = n*(in_levels+1) + cnt101
    ELSE
      totalnumberofelements = n*(in_levels+3) + cnt101
    END IF

    IF (PreserveBaseline) &
        totalnumberofelements = totalnumberofelements + Mesh_in % NumberOfBoundaryElements
    ALLOCATE(Mesh_out % Elements(totalnumberofelements))
    
    ! Generate volume bulk elements:
    ! ------------------------------

    Mesh_out % MaxElementNodes = 0

    NeedEdges=.FALSE.
    n=Mesh_in % NumberOfNodes
    cnt=0; dg_n  = 0
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBulkElements

        cnt=cnt+1
        Mesh_out % Elements(cnt) = Mesh_in % Elements(j)

        l_n=0
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+i*n
        END DO
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          IF( Rotate2Pi .AND. i==in_levels ) THEN
            ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)
          ELSE
            ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+(i+1)*n
          END IF
        END DO
        Mesh_out % Elements(cnt) % NDOFs = l_n
        Mesh_out % MaxElementNodes=MAX(Mesh_out % MaxElementNodes,l_n)

        SELECT CASE(l_n)
        CASE(6)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(706)
        CASE(8)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(808)
        END SELECT

        Mesh_out % Elements(cnt) % GElementIndex = &
             Mesh_in % Elements(j) % GelementIndex + gelements*i

        Mesh_out % Elements(cnt) % ElementIndex = cnt
        ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n)) 
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % NodeIndexes = ind(1:l_n)
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    END DO
    Mesh_out % NumberOfBulkElements=cnt

    max_bid=0
    max_baseline_bid=0

    ! include edges (see below)
    NeedEdges =  (NeedEdges .OR. PreserveEdges)
    
    ! -------------------------------------------------------
    IF (PreserveBaseline) THEN
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

        ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
        Mesh_out % Elements(cnt) % BoundaryInfo = &
           Mesh_in % Elements(k) % BoundaryInfo

        max_bid = MAX(max_bid, Mesh_in % Elements(k) % &
                BoundaryInfo % Constraint)

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in %  NumberOfBulkElements*(in_levels+1)+ &
	                   (in_levels+2)*Mesh_in % NumberOfBoundaryElements+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Right => &
              Mesh_out % Elements(Mesh_in % NumberOfBulkElements*(in_levels+1)+ &
	      (in_levels+2)*Mesh_in % NumberOfBoundaryElements+l)
        END IF

        IF(Mesh_in % Elements(k) % TYPE % ElementCode>=200) THEN
          Mesh_out % Elements(cnt) % NDOFs = 2
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(2)) 
          ind(1) = Mesh_in % Elements(k) % NodeIndexes(1)
          ind(2) = Mesh_in % Elements(k) % NodeIndexes(2)
          Mesh_out % Elements(cnt) % NodeIndexes = ind(1:2)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(202)
        ELSE
          Mesh_out % Elements(cnt) % NDOFs = 1
          l=SIZE(Mesh_in % Elements(k) % NodeIndexes)
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l))
          Mesh_out % Elements(cnt) % NodeIndexes = &
            Mesh_in % Elements(k) % NodeIndexes
          Mesh_out % Elements(cnt) % TYPE => &
             Mesh_in % Elements(k) % TYPE
        END IF
        Mesh_out % Elements(cnt) % DGDOFs = 0
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % ElementIndex = cnt
        Mesh_out % Elements(cnt) % PDefs => NULL()
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    
      IF(isParallel) THEN
        j=max_bid
        CALL MPI_ALLREDUCE(j,max_bid,1, &
            MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
      END IF

      max_baseline_bid = max_bid

    END IF


    ! Add side boundaries with the bottom mesh boundary id's:
    ! (or shift ids if preserving the baseline boundary)
    ! -------------------------------------------------------
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

        ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
        Mesh_out % Elements(cnt) % BoundaryInfo = &
           Mesh_in % Elements(k) % BoundaryInfo

        Mesh_out % Elements(cnt) % BoundaryInfo % constraint = &
           Mesh_out % Elements(cnt) % BoundaryInfo % constraint + max_baseline_bid

        max_bid = MAX(max_bid, max_baseline_bid + &
           Mesh_in % Elements(k) % BoundaryInfo % Constraint)

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        IF(Mesh_in % Elements(k) % TYPE % ElementCode>=200) THEN
          Mesh_out % Elements(cnt) % NDOFs = 4
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(4)) 

          ind(1) = Mesh_in % Elements(k) % NodeIndexes(1)+i*n
          ind(2) = Mesh_in % Elements(k) % NodeIndexes(2)+i*n

          IF( Rotate2Pi .AND. i==in_levels ) THEN
            ind(3) = Mesh_in % Elements(k) % NodeIndexes(2)
            ind(4) = Mesh_in % Elements(k) % NodeIndexes(1)
          ELSE
            ind(3) = Mesh_in % Elements(k) % NodeIndexes(2)+(i+1)*n
            ind(4) = Mesh_in % Elements(k) % NodeIndexes(1)+(i+1)*n
          END IF
            Mesh_out % Elements(cnt) % NodeIndexes = ind(1:4)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(404)
        ELSE
          Mesh_out % Elements(cnt) % NDOFs = 1
          l=SIZE(Mesh_in % Elements(k) % NodeIndexes)
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l))
          Mesh_out % Elements(cnt) % NodeIndexes = &
            Mesh_in % Elements(k) % NodeIndexes+i*n
          Mesh_out % Elements(cnt) % TYPE => &
             Mesh_in % Elements(k) % TYPE
        END IF 
        Mesh_out % Elements(cnt) % ElementIndex = cnt
        Mesh_out % Elements(cnt) % DGDOFs = 0
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % PDefs => NULL()
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    END DO

    !Take care of extra 101 elements
    !-------------------------------

    IF(cnt101 > 0) THEN
       DO j=1,Mesh_in % NumberOfBoundaryElements
          k = j + Mesh_in % NumberOfBulkElements

          IF(Mesh_in % Elements(k) % TYPE % ElementCode /= 101) CYCLE
          cnt=cnt+1

          Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

          ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
          Mesh_out % Elements(cnt) % BoundaryInfo = &
               Mesh_in % Elements(k) % BoundaryInfo

          Mesh_out % Elements(cnt) % BoundaryInfo % constraint = &
               Mesh_out % Elements(cnt) % BoundaryInfo % constraint + max_baseline_bid

          max_bid = MAX(max_bid, max_baseline_bid + &
               Mesh_in % Elements(k) % BoundaryInfo % Constraint)

          Mesh_out % Elements(cnt) % NDOFs = 1
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(1))
          Mesh_out % Elements(cnt) % NodeIndexes = &
               Mesh_in % Elements(k) % NodeIndexes+(in_levels+1)*n
          Mesh_out % Elements(cnt) % TYPE => &
               Mesh_in % Elements(k) % TYPE

          Mesh_out % Elements(cnt) % ElementIndex = cnt
          Mesh_out % Elements(cnt) % DGDOFs = 0
          Mesh_out % Elements(cnt) % DGIndexes => NULL()
          Mesh_out % Elements(cnt) % PDefs => NULL()
          Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
          Mesh_out % Elements(cnt) % FaceIndexes => NULL()
          Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
       END DO
    END IF
    
    IF(isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    WRITE( Message,'(A,I0)') 'First Extruded BC set to: ',max_bid+1
    CALL Info('MeshExtrude',Message,Level=8)

    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    WRITE( Message,'(A,I0)') 'Number of new BCs for layers: ',max_body
    CALL Info('MeshExtrude',Message,Level=8)


    ! Add start and finish planes except if we have a full rotational symmetry
    IF( .NOT. Rotate2Pi ) THEN

    ! Add bottom boundary:
    ! --------------------
    DO i=1,Mesh_in % NumberOfBulkElements
      cnt=cnt+1

      Mesh_out % Elements(cnt) = Mesh_in % Elements(i)

      l_n=Mesh_in % Elements(i) % TYPE % NumberOfNodes
      Mesh_out % Elements(cnt) % NDOFs = l_n

      ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
      Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
           Mesh_out % Elements(i)
      Mesh_out % Elements(cnt) % BoundaryInfo % Right => NULL()

      bcid = max_bid + Mesh_out % Elements(cnt) % BodyId
      Mesh_out % Elements(cnt) % BoundaryInfo % Constraint = bcid

      Mesh_out % Elements(cnt) % BodyId = 0
      IF( bcid<=CurrentModel % NumberOfBCs) THEN
        j=ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
        IF(Found) Mesh_out % Elements(cnt) % BodyId=j
      END IF

      ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n))
      Mesh_out % Elements(cnt) % NodeIndexes = &
        Mesh_in % Elements(i) % NodeIndexes
      Mesh_out % Elements(cnt) % ElementIndex = cnt
      Mesh_out % Elements(cnt) % TYPE => &
        Mesh_in % Elements(i) % TYPE
      Mesh_out % Elements(cnt) % DGDOFs = 0
      Mesh_out % Elements(cnt) % DGIndexes => NULL()
      Mesh_out % Elements(cnt) % PDefs => NULL()
      Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
      Mesh_out % Elements(cnt) % FaceIndexes => NULL()
      Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
    END DO

    ! Add top boundary:
    ! -----------------
    DO i=1,Mesh_in % NumberOfBulkElements
      cnt=cnt+1

      Mesh_out % Elements(cnt) = Mesh_in % Elements(i)

      l_n=Mesh_in % Elements(i) % TYPE % NumberOfNodes
      Mesh_out % Elements(cnt) % NDOFs = l_n

      ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
      Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
           Mesh_out % Elements(in_levels*Mesh_in % NumberOfBulkElements+i)
      Mesh_out % Elements(cnt) % BoundaryInfo % Right => NULL()

      bcid = max_bid + Mesh_out % Elements(cnt) % BodyId + max_body
      Mesh_out % Elements(cnt) % BoundaryInfo % Constraint = bcid

      Mesh_out % Elements(cnt) % BodyId = 0
      IF( bcid<=CurrentModel % NumberOfBCs) THEN
        j=ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
        IF(Found) Mesh_out % Elements(cnt) % BodyId=j
      END IF

      ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n))
      Mesh_out % Elements(cnt) % NodeIndexes = &
        Mesh_in % Elements(i) % NodeIndexes+(in_Levels+1)*n
      Mesh_out % Elements(cnt) % ElementIndex = cnt
      Mesh_out % Elements(cnt) % TYPE => &
        Mesh_in % Elements(i) % TYPE
      Mesh_out % Elements(cnt) % DGDOFs = 0
      Mesh_out % Elements(cnt) % DGIndexes => NULL()
      Mesh_out % Elements(cnt) % PDefs => NULL()
      Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
      Mesh_out % Elements(cnt) % FaceIndexes => NULL()
      Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
    END DO

    END IF ! .NOT. Rotate2Pi
    

    Mesh_out % NumberOfBoundaryElements=cnt-Mesh_out % NumberOfBulkElements

    Mesh_out % Name=Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs  = Mesh_out % MaxElementNodes
    Mesh_out % Stabilize = Mesh_in % Stabilize
    Mesh_out % MeshDim = 3
    CurrentModel % Dimension = 3

    CALL PrepareMesh( CurrentModel, Mesh_out, isParallel )
    
    ExtrudedMeshName = ListGetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
    IF(Found) THEN
      CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
    END IF

!------------------------------------------------------------------------------
  END FUNCTION MeshExtrude
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> As the previous one except the extrusion is done in parallel for single meshes
!> that each take an internal in the extruded direction. This affects the coordinates
!> but also the communication pattern. A separate routine was made in order to avoid
!> introducing of bugs as the internal extrusion is a widely used feature. 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrudeSlices(Mesh_in, in_levels) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh_in, Mesh_out
    INTEGER :: in_levels
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: ExtrudedMeshName
    INTEGER :: i,j,k,l,n,m,cnt,ind(8),bid,max_bid,l_n,max_body,bcid,&
        ExtrudedCoord,dg_n,totalnumberofelements
    INTEGER, POINTER :: pInds(:)
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    TYPE(Element_t), POINTER :: Element
    INTEGER :: nnodes,gnodes,gelements,ierr,nlev,ilev,&
        nParMesh,nParExt,OrigPart,ElemCode,bodyid
    LOGICAL :: isParallel, SingleIn, Found, TopBC, BotBC
    INTEGER,ALLOCATABLE :: ChildBCs(:)
    REAL(KIND=dp)::w,MinCoord,MaxCoord,CurrCoord,zmin,zmax
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
!------------------------------------------------------------------------------

    ! The historical choice in_levels in annoying when we want to split the divisions.
    nlev = in_levels+1
    
    CALL Info('MeshExtrudeSlices','Creating '//TRIM(I2S(nlev))//' extruded element layers',Level=10)

    IF( ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found ) ) &
        CALL Fatal('MeshExtrudeSlices','The slice version cannot handle "Preserve Baseline"!')
    
    IF( ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Rotational',Found ) ) &
        CALL Fatal('MeshExtrudeSlices','The slice version cannot handle "Extruded Mesh Rotational"!')    
    
    isParallel = ( ParEnv % PEs > 1 )
    SingleIn = Mesh_in % SingleMesh
    
    ! Create the division for the 1D unit mesh
    !--------------------------------------------
    ALLOCATE( Wtable( 0: nlev ) )
    CALL UnitSegmentDivision( Wtable, nlev )
    
    ! In parallel let us pick only our own share of the
    ! division. This logic makes it possible to have nonuniform divisions easily.
    ! The number of element layers is evenly distributed among partitions. 
    !-------------------------------------------------------------------------------
    IF( isParallel ) THEN
      nParExt = ParEnv % PEs 
      nParMesh = ListGetInteger( CurrentModel % Simulation,'Parallel Mesh Modulo',Found)
      IF(.NOT. Found) THEN
        nParMesh = 1
        IF(.NOT. SingleIn ) THEN
          CALL Fatal('MeshExtrudedSlices','This routine expects either Mesh Modulo or Single Mesh!')
        END IF
      END IF
      
      nParExt = nParExt / nParMesh                    
      IF( MODULO(nlev,nParExt) /= 0 ) THEN
        CALL Fatal('MeshExtrudedSlices','Number of element layers '//TRIM(I2S(nlev))//&
            ' not divisible by '//TRIM(I2S(ParEnv % PEs)))
      END IF
      nlev = nlev / nParExt
      IF(nlev < 2) THEN
        CALL Fatal('MeshExtrudedSlices','At least two element layers needed in each partition!')
      END IF
      ilev = ( ParEnv % MyPe / nParMesh ) * nlev
      Wtable(0:nlev) = Wtable(ilev:nlev+ilev) 
    ELSE
      nParExt = 1
      nParMesh = 1 
    END IF
        
    ! Allocate extruded mesh:
    ! We do this only after splitting the division.
    ! ---------------------------------------------
    n = Mesh_in % NumberOfNodes
    nnodes = (nlev+1)*n

    Mesh_out => AllocateMesh()
    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )
    
    gnodes = Mesh_in % NumberOfNodes
    gelements = Mesh_in % NumberOfBulkElements

    Mesh_out % SingleMesh = .FALSE.
    
    ExtrudedCoord = ListGetInteger( CurrentModel % Simulation,'Extruded Coordinate Index', &
        Found, minv=1,maxv=3 )
    IF(.NOT. Found) ExtrudedCoord = 3 
    
    IF( ExtrudedCoord == 1 ) THEN
      ActiveCoord => Mesh_out % Nodes % x
    ELSE IF( ExtrudedCoord == 2 ) THEN
      ActiveCoord => Mesh_out % Nodes % y
    ELSE IF( ExtrudedCoord == 3 ) THEN
      ActiveCoord => Mesh_out % Nodes % z
    END IF

    MinCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Min Coordinate',Found )
    IF(.NOT. Found) MinCoord = 0.0_dp
    MaxCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Max Coordinate',Found )
    IF(.NOT. Found) MaxCoord = 1.0_dp

     
    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo

      IF(.NOT. ASSOCIATED( PI_in ) ) CALL Fatal('MeshExtrudeSlices','PI_in not associated!')
      IF(.NOT. ASSOCIATED( PI_out ) ) CALL Fatal('MeshExtrudeSlices','PI_out not associated!')
            
      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % NodeInterface(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      IF(.NOT. SingleIn ) THEN
        IF(.NOT. ASSOCIATED( PI_in % NeighbourList ) ) THEN
          CALL Fatal('MeshExtrudeSlices','Neighnours not associated!')
        END IF
      END IF
        
      IF(.NOT. SingleIn ) THEN
        ! Count own nodes
        j=0
        DO i=1,Mesh_in % NumberOfNodes
          IF(.NOT. ASSOCIATED(PI_in % NeighbourList(i) % Neighbours ) ) THEN
            j = j + 1
          ELSE IF (PI_in % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
            j=j+1
          END IF
        END DO
        CALL MPI_ALLREDUCE(j,gnodes,1, &
            MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
        gnodes = gnodes / nParExt
        
        j=0
        DO i=1,Mesh_in % NumberOfBulkElements
          IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
        END DO
        CALL MPI_ALLREDUCE(j,gelements,1, &
            MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
        gelements = gelements / nParExt

        !PRINT *,'nParExt:',ParEnv % Mype, nParExt, nParMesh, gnodes,gelements        
      END IF
    END IF

    CALL Info('MeshExtrudeSlices','Number of nodes in layer: '//TRIM(I2S(gnodes)),Level=12)
    CALL Info('MeshExtrudeSlices','Number of elements in layer: '//TRIM(I2S(gelements)),Level=12)
    
    !CALL Info('MeshExtrudeSlices','Number of extruded nodes: '//TRIM(I2S((nlev+1)*gnodes)),Level=7)
    !CALL Info('MeshExtrudeSlices','Number of exruded elements: '//TRIM(I2S(nlev*gelements)),Level=7)
    
    cnt=0
    DO i=0,nlev

      w = Wtable( i ) 
      CurrCoord = w * MaxCoord + (1-w) * MinCoord      
      
      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          m = 1
          IF( nParMesh > 1 ) THEN
            IF( ASSOCIATED( PI_in % NeighbourList(j) % Neighbours ) ) THEN
              m = SIZE(PI_in % NeighbourList(j) % Neighbours)
            END IF
          END IF
          IF(i==0 .AND. ParEnv % MyPe > (nParMesh-1) ) THEN            
            k = 2*m
          ELSE IF(i==nlev .AND. ParEnv % MyPe < ParEnv % PEs- nParMesh ) THEN
            k = 2*m
          ELSE
            k = m
          END IF

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(k))
          PI_out % NodeInterface(cnt) = (k>1)
        
          DO k=1,m
            IF(m>1) THEN
              OrigPart = PI_in % NeighbourList(j) % Neighbours(k)
            ELSE
              OrigPart = ParEnv % MyPe
            END IF                       

            IF(SingleIn) THEN
              l = j + (ilev+i) * gnodes
            ELSE
              l = MODULO(PI_in % GlobalDOFs(j)-1,gnodes)+1 + (ilev+i) * gnodes 
            END IF
            PI_out % GlobalDOFs(cnt) = l
                                     
            IF(i==0 .AND. ParEnv % MyPe > nParMesh-1 ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(2*k-1) = OrigPart
              PI_out % NeighbourList(cnt) % Neighbours(2*k) = OrigPart-1            
            ELSE IF(i==nlev .AND. ParEnv % MyPe < ParEnv % PEs-nParMesh ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(2*k-1) = OrigPart+1
              PI_out % NeighbourList(cnt) % Neighbours(2*k) = OrigPart
            ELSE
              PI_out % NeighbourList(cnt) % Neighbours(k) = OrigPart 
            END IF                       
          END DO
          
        END IF
      END DO
    END DO
    
    Mesh_out % NumberOfNodes = cnt
    Mesh_out % Nodes % NumberOfNodes = cnt

    ! Calculate exactly and allocate the number of extruded elements
    n = Mesh_in % NumberOfBulkElements + Mesh_in % NumberOfBoundaryElements
    totalnumberofelements = n*nlev
    IF( ParEnv % MyPe < nParMesh ) totalnumberofelements = &
        totalnumberofelements + Mesh_in % NumberOfBulkElements 
    IF( ParEnv % MyPe >= ParEnv % PEs-nParMesh ) totalnumberofelements = &
        totalnumberofelements + Mesh_in % NumberOfBulkElements 
      
    ALLOCATE(Mesh_out % Elements(totalnumberofelements))

    
    ! Generate volume bulk elements:
    ! ------------------------------
    Mesh_out % MaxElementNodes = 0
    n = Mesh_in % NumberOfNodes
    cnt=0; dg_n  = 0

    DO i=0,nlev-1
      DO j=1,Mesh_in % NumberOfBulkElements

        cnt = cnt+1
        Element => Mesh_out % Elements(cnt)
        Element = Mesh_in % Elements(j)

        l_n=0
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+i*n
        END DO
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+(i+1)*n
        END DO
        Element % NDOFs = l_n
        Mesh_out % MaxElementNodes = MAX(Mesh_out % MaxElementNodes,l_n)

        SELECT CASE(l_n)
        CASE(6)
          Element % TYPE => GetElementType(706)
        CASE(8)
          Element % TYPE => GetElementType(808)
        END SELECT

        IF( isParallel ) THEN
          IF(SingleIn) THEN
            l = j + (ilev+i) * gelements
          ELSE
            l = MODULO(Mesh_in % Elements(j) % GElementIndex-1,gelements)+1 + (ilev+i) * gelements 
          END IF
          Element % GElementIndex = l
        END IF
          
        Element % ElementIndex = cnt
        ALLOCATE(Element % NodeIndexes(l_n)) 
        Element % NodeIndexes = ind(1:l_n)
      END DO
    END DO
    Mesh_out % NumberOfBulkElements = cnt

    
    ! Add side boundaries with the bottom mesh boundary id's:
    ! -------------------------------------------------------
    max_bid = 0
    DO i=0,nlev-1
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Element => Mesh_out % Elements(cnt)
        Element = Mesh_in % Elements(k)        
        ALLOCATE(Element % BoundaryInfo)

        Element % BoundaryInfo = Mesh_in % Elements(k) % BoundaryInfo

        bid = Mesh_in % Elements(k) % BoundaryInfo % Constraint
        max_bid = MAX(max_bid, bid )

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l = Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Element % BoundaryInfo % Left => &
              Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l = Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Element % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        ElemCode = Mesh_in % Elements(k) % TYPE % ElementCode        
        m = 2*MODULO(ElemCode,100)        
        Element % NDOFs = m
        ALLOCATE(Element % NodeIndexes(m))
        pInds => Element % NodeIndexes
               
        IF(ElemCode == 202) THEN
          pInds(1) = Mesh_in % Elements(k) % NodeIndexes(1)+i*n
          pInds(2) = Mesh_in % Elements(k) % NodeIndexes(2)+i*n
          pInds(3) = Mesh_in % Elements(k) % NodeIndexes(2)+(i+1)*n
          pInds(4) = Mesh_in % Elements(k) % NodeIndexes(1)+(i+1)*n
          Mesh_out % Elements(cnt) % TYPE => GetElementType(404)
        ELSE IF(ElemCode == 101 ) THEN
          pInds(1) = Mesh_in % Elements(k) % NodeIndexes(1) +i*n
          pInds(2) = Mesh_in % Elements(k) % NodeIndexes(1) +(i+1)*n
        ELSE
          CALL Fatal('MeshExtrudeSlices','Cannot extrude boundary element: '//TRIM(I2S(ElemCode)))
        END IF
        Mesh_out % Elements(cnt) % ElementIndex = cnt
      END DO
    END DO

    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    
    WRITE( Message,'(A,I0)') 'First Extruded BC set to: ',max_bid+1
    CALL Info('MeshExtrudeSlices',Message,Level=8)

    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    WRITE( Message,'(A,I0)') 'Number of new BCs for each layer: ',max_body
    CALL Info('MeshExtrudeSlices',Message,Level=8)

    ALLOCATE(ChildBCs(2*max_body))
    ChildBCs = -1
           
    ! Add bottom boundary:
    ! --------------------
    IF( ParEnv % PEs == 1 .OR. ParEnv % MyPe < nParMesh ) THEN  
      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        Element => Mesh_out % Elements(cnt) 
        
        Element = Mesh_in % Elements(i)

        l_n = Mesh_in % Elements(i) % TYPE % NumberOfNodes
        Element % NDOFs = l_n

        ALLOCATE(Element % BoundaryInfo)
        Element % BoundaryInfo % Left => Mesh_out % Elements(i)
        Element % BoundaryInfo % Right => NULL()

        bodyid = Mesh_in % Elements(i) % BodyId                
        bcid = max_bid + bodyid
        Element % BoundaryInfo % Constraint = bcid

        ChildBCs(2*bodyid-1) = bcid 

        Element % BodyId = 0
        IF( bcid <= CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
          IF(Found) Element % BodyId = j
        END IF

        ALLOCATE(Element % NodeIndexes(l_n))
        Element % NodeIndexes = Mesh_in % Elements(i) % NodeIndexes
        Element % ElementIndex = cnt
        Element % TYPE => Mesh_in % Elements(i) % TYPE
      END DO
    END IF

    
    ! Add top boundary:
    ! -----------------
    IF( ParEnv % PEs == 1 .OR. ParEnv % MyPe >= ParEnv % PEs - nParMesh ) THEN
      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        Element => Mesh_out % Elements(cnt) 
        
        Element = Mesh_in % Elements(i)

        l_n = Mesh_in % Elements(i) % TYPE % NumberOfNodes
        Element % NDOFs = l_n

        ALLOCATE(Element % BoundaryInfo)
        Element % BoundaryInfo % Left => &
            Mesh_out % Elements((nlev-1)*Mesh_in % NumberOfBulkElements+i)
        Element % BoundaryInfo % Right => NULL()
        
        bodyid = Mesh_in % Elements(i) % BodyId                
        bcid = max_bid + bodyid + max_body
        Element % BoundaryInfo % Constraint = bcid

        ChildBCs(2*bodyid) = bcid 
        
        Element % BodyId = 0
        IF( bcid<=CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
          IF(Found) Element % BodyId = j
        END IF

        ALLOCATE(Element % NodeIndexes(l_n))
        Element % NodeIndexes = Mesh_in % Elements(i) % NodeIndexes+nlev*n
        Element % ElementIndex = cnt
        Element % TYPE => Mesh_in % Elements(i) % TYPE
      END DO
    END IF
    
    IF( cnt /= totalnumberofelements ) THEN
      CALL Fatal('MeshExtrudedSlices','Mismatch between allocated and set elements: '//&
          TRIM(I2S(totalnumberofelements))//' vs. '//TRIM(I2S(cnt)))
    END IF

    ! Set some unset stuff to be on the safe side
    DO i=1,cnt
      Element => Mesh_out % Elements(i)
      Element % DGDOFs = 0
      Element % DGIndexes => NULL()
      Element % PDefs => NULL()
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()
      Element % BubbleIndexes => NULL()
    END DO
         
    Mesh_out % NumberOfBoundaryElements = cnt - Mesh_out % NumberOfBulkElements
    
    Mesh_out % Name = Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs = Mesh_out % MaxElementNodes
    Mesh_out % Stabilize = Mesh_in % Stabilize
    Mesh_out % MeshDim = 3
    CurrentModel % DIMENSION = 3


    ! Let us mark the child BCs to the bodies that they originate from.
    BLOCK
      INTEGER, POINTER :: TmpPair(:), TmpBCs(:) 
      TYPE(ValueList_t), POINTER :: vList

      ALLOCATE(TmpBCs(2*max_body))
      TmpBCs = ChildBCs

      IF( ParEnv % PEs > 1 ) THEN
        CALL MPI_ALLREDUCE(TmpBCs,ChildBCs,2*max_body, &
            MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
      END IF

      DO i=1,CurrentModel % NumberOfBodies
        vList => CurrentModel % Bodies(i) % Values
        IF( ASSOCIATED(vList) ) THEN
          NULLIFY(TmpPair)
          ALLOCATE(TmpPair(2))
          TmpPair(1) = ChildBCs(2*i-1)
          TmpPair(2) = ChildBCs(2*i)
          CALL ListAddIntegerArray(vList,'Extruded Child BCs',2,TmpPair)

          IF( InfoActive(20) ) THEN
            PRINT *,'Extruded Child BCs for body:',i,TmpPair
          END IF
          NULLIFY(TmpPair)
        END IF
      END DO

      DEALLOCATE(TmpBCs)
    END BLOCK
      
    
    ExtrudedMeshName = ListGetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
    IF(Found) THEN
      IF( ParEnv % PEs == 1 ) THEN
        CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
      ELSE
        CALL WriteMeshToDisk2(CurrentModel, Mesh_out, ExtrudedMeshName, ParEnv % MyPe )
      END IF
    END IF

    CALL PrepareMesh( CurrentModel, Mesh_out, isParallel )
    
!------------------------------------------------------------------------------
  END FUNCTION MeshExtrudeSlices
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Writes the mesh to disk. Note that this does not include the information
!> of shared nodes needed in parallel computation. This may be used for 
!> debugging purposes and for adaptive solution, for example. 
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk( NewMesh, Path )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: Path
    TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,MaxNodes,ElmCode,Parent1,Parent2
!------------------------------------------------------------------------------

    OPEN( 1,FILE=TRIM(Path) // '/mesh.header',STATUS='UNKNOWN' )
    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, NewMesh % NumberOfBoundaryElements
    
    WRITE( 1,'(i0)' ) 2
    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(k) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(k) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(k) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBoundaryElements

    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(i) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(i) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(i) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBulkElements
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.nodes', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       WRITE(1,'(i0,a,3e23.15)',ADVANCE='NO') i,' -1 ', &
            NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.elements', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       WRITE(1,'(3(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % Elements(i) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.boundary', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex
       WRITE(1,'(5(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(k) % BoundaryInfo % Constraint, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % Elements(k) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)
!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk2(Model, NewMesh, Path, Partition )
!------------------------------------------------------------------------------
    USE Types
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: NewMesh
    CHARACTER(LEN=*) :: Path
    INTEGER, OPTIONAL :: Partition
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeList(100),ElmCodeCounts(100),&
         Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared
    INTEGER, POINTER :: BList(:)
    INTEGER, ALLOCATABLE :: ElementCodes(:)
    LOGICAL :: Parallel, WarnNoTarget, Found
    CHARACTER(LEN=MAX_NAME_LEN) :: headerFN, elementFN, nodeFN,&
         boundFN, sharedFN
!------------------------------------------------------------------------------

    IF(PRESENT(Partition)) THEN
       Parallel = .TRUE.
       WRITE(headerFN, '(A,I0,A)') '/part.',Partition+1,'.header'
       WRITE(elementFN, '(A,I0,A)') '/part.',Partition+1,'.elements'
       WRITE(nodeFN, '(A,I0,A)') '/part.',Partition+1,'.nodes'
       WRITE(boundFN, '(A,I0,A)') '/part.',Partition+1,'.boundary'
       WRITE(sharedFN, '(A,I0,A)') '/part.',Partition+1,'.shared'
    ELSE
       Parallel = .FALSE.
       headerFN = '/mesh.header'
       elementFN = '/mesh.elements'
       nodeFN = '/mesh.nodes'
       boundFN = '/mesh.boundary'
    END IF

    !Info for header file

    ElmCodeList = 0 !init array
    NumElmCodes = 0
    NumElements = NewMesh % NumberOfBoundaryElements + &
         NewMesh % NumberOfBulkElements
    ALLOCATE(ElementCodes(NumElements))

    !cycle to bring element code list into array-inquirable form
    DO i=1,NumElements
       ElementCodes(i) = NewMesh % Elements(i) % TYPE % ElementCode
    END DO

    DO i=NumElements,1,-1 !this should give element codes increasing value, which appears to be
                          !'standard' though I doubt it matters
       IF(ANY(ElmCodeList == ElementCodes(i))) CYCLE
       NumElmCodes = NumElmCodes + 1
       ElmCodeList(NumElmCodes) = ElementCodes(i)
    END DO

    DO j=1,NumElmCodes
       ElmCodeCounts(j) = COUNT(ElementCodes == ElmCodeList(j))
    END DO

    !Write header file
    OPEN( 1,FILE=TRIM(Path) // headerFN,STATUS='UNKNOWN' )
    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, &
         NewMesh % NumberOfBoundaryElements

    WRITE( 1,'(i0)' ) NumElmCodes
    DO j=1,NumElmCodes
       WRITE( 1,'(i0,x,i0,x)' ) ElmCodeList(j),ElmCodeCounts(j)
    END DO
    IF(Parallel) THEN !need number of shared nodes
       NoShared = 0
       DO i=1,NewMesh % NumberOfNodes
          IF(SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
               Neighbours) > 1) THEN
             NoShared = NoShared + 1
          END IF
       END DO
       WRITE( 1,'(i0,x,i0)') NoShared, 0
    END IF
    CLOSE(1)

    !Write nodes file
    OPEN( 1,FILE=TRIM(Path) // nodeFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       IF (Parallel) THEN
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % ParallelInfo % GlobalDOFs(i)
       ELSE
          WRITE(1,'(i0,x)', ADVANCE='NO') i
       END IF
       WRITE(1,'(a,x,ES17.10,x,ES17.10,x,ES17.10)',ADVANCE='NO') &
            ' -1 ', NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    !Write elements file
    OPEN( 1,FILE=TRIM(Path) // elementFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       IF(Parallel) THEN
          ElemID = NewMesh % Elements(i) % GElementIndex
       ELSE
          ElemID = i
       END IF
       WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') ElemID, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(i) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(i) % NodeIndexes(j)
          END IF
          WRITE(1,'(i0,x)', ADVANCE='NO') m
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    !Write boundary file
    WarnNoTarget = .FALSE.
    OPEN( 1,FILE=TRIM(Path) // boundFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex

       IF(Parallel) THEN
          IF(parent1 /= 0) parent1 = NewMesh % Elements(parent1) % GElementIndex
          IF(parent2 /= 0) parent2 = NewMesh % Elements(parent2) % GElementIndex
       END IF

       Constraint = NewMesh % Elements(k) % BoundaryInfo % Constraint
       BList => ListGetIntegerArray( Model % BCs(Constraint) % Values, &
            'Target Boundaries', Found )
       IF(Found) THEN
          IF(SIZE(BList) > 1) THEN
             CALL WARN("WriteMeshToDisk2",&
                  "A BC has more than one Target Boundary, SaveMesh output will not match input!")
          END IF
          meshBC = BList(1)
       ELSE
          WarnNoTarget = .TRUE.
          meshBC = Constraint
       END IF

       !This meshBC stuff will *only* work if each BC has only 1 target boundary
       WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            meshBC, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(k) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(k) % NodeIndexes(j)
          END IF
          WRITE(1,'(x,i0)', ADVANCE='NO') m
       END DO
       WRITE(1,*) !blank write statement to create new line without extra space.
    END DO
    CLOSE(1)

    IF(WarnNoTarget) THEN
       CALL WARN("WriteMeshToDisk2","Couldn't find a Target Boundary, assuming mapping to self")
    END IF

    IF(.NOT. Parallel) RETURN

    !Write .shared file
    !Need to create part.n.shared from Mesh % ParallelInfo %
    !NeighbourList % Neighbours.
    OPEN( 1,FILE=TRIM(Path) // sharedFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       nneigh = SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
            Neighbours)
       IF(nneigh < 2) CYCLE
       WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') &
            NewMesh % ParallelInfo % GlobalDOFs(i),nneigh
       DO j=1,nneigh
          WRITE(1,'(I0, x)',ADVANCE='NO') NewMesh % ParallelInfo %&
               NeighbourList(i) % Neighbours(j) + 1
       END DO
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)


!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk2
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDiskPartitioned(Model, Mesh, Path, &
      ElementPart, NeighbourList )
!------------------------------------------------------------------------------
    USE Types
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: Path
    INTEGER, POINTER :: ElementPart(:)
    TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: NoBoundaryElements, NoBulkElements, NoNodes, NoPartitions, Partition
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeCounts(827),&
         Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared
    LOGICAL :: Found, Hit
    CHARACTER(LEN=MAX_NAME_LEN) :: DirectoryName, PrefixName
!------------------------------------------------------------------------------

    NoPartitions = MAXVAL( ElementPart ) 
    NumElmCodes = 0
    NumElements = Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements
        
    WRITE(DirectoryName, '(A,A,I0)') TRIM(PATH),'/partitioning.',NoPartitions
    CALL MakeDirectory( TRIM(DirectoryName) // CHAR(0) )
    CALL Info('WriteMeshToDiskPartitioned','Writing parallel mesh to disk: '//TRIM(DirectoryName))
   

    DO Partition = 1, NoPartitions 
      
      CALL Info('WriteMeshToDiskPartitioned','Writing piece to file: '//TRIM(I2S(Partition)),Level=12)
      
      WRITE( PrefixName,'(A,A,I0)') TRIM(DirectoryName),'/part.',Partition  

      CALL Info('WriteMeshToDiskPartitioned','Write nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.nodes', STATUS='UNKNOWN' )
      NoNodes = 0
      DO i=1,Mesh % NumberOfNodes
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          WRITE(1,'(I0,x,I0,x,3ES17.10)') i,-1, &
              Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)
          NoNodes = NoNodes + 1
        END IF
      END DO
      CLOSE(1)
      

      CALL Info('WriteMeshToDiskPartitioned','Write shared nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.shared', STATUS='UNKNOWN' )
      NoShared = 0
      DO i=1,Mesh % NumberOfNodes
        nneigh = SIZE( NeighbourList(i) % Neighbours )
        IF( nneigh <= 1 ) CYCLE
        
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          NoShared = NoShared + 1
          WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') i,nneigh
          DO j=1,nneigh
            WRITE(1,'(I0, x)',ADVANCE='NO') NeighbourList(i) % Neighbours(j) 
          END DO
          WRITE( 1,* ) ''
        END IF
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write elements file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.elements', STATUS='UNKNOWN' )
      NoBulkElements = 0
      ElmCodeCounts = 0      
      DO i=1,Mesh % NumberOfBulkElements
        IF( ElementPart(i) /= Partition ) CYCLE

        Element => Mesh % Elements(i)
        WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') i, &
            Element % BodyId, Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) ''
        
        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBulkElements = NoBulkElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write boundary file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.boundary', STATUS='UNKNOWN' )
      NoBoundaryElements = 0
      DO i=Mesh % NumberOfBulkElements +1 ,&
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(i)
       
        parent1 = 0
        parent2 = 0
        Constraint = 0
        
        IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
          IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) &
              parent1 = Element % BoundaryInfo % Left % ElementIndex
          IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) &
              parent2 = Element % BoundaryInfo % Right % ElementIndex        
          Constraint = Element % BoundaryInfo % Constraint
        END IF

        Hit = .FALSE.
        IF( parent1 > 0 ) THEN
          IF( ElementPart( parent1 ) == Partition ) Hit = .TRUE.
        END IF
        IF( parent2 > 0 ) THEN
          IF( ElementPart( parent2 ) == Partition ) Hit = .TRUE.
        END IF

        IF( .NOT. Hit ) CYCLE

        WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            Constraint, Parent1, Parent2,&
            Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(x,i0)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) 

        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBoundaryElements = NoBoundaryElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write header file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.header',STATUS='UNKNOWN' )
      NumElmCodes = COUNT( ElmCodeCounts > 0 ) 
      WRITE( 1,'(i0,x,i0,x,i0)' ) NoNodes, &
          NoBulkElements, NoBoundaryElements      
      WRITE( 1,'(i0)' ) NumElmCodes
      DO i=SIZE(ElmCodeCounts),1,-1
        IF( ElmCodeCounts(i) == 0 ) CYCLE
        WRITE( 1,'(i0,x,i0,x)' ) i,ElmCodeCounts(i)
      END DO
      WRITE( 1,'(i0,x,i0)') NoShared, 0
      CLOSE(1)
      
      CALL Info('WriteMeshToDiskPartitioned','Done writing partition',Level=12)
    END DO

    CALL Info('WriteMeshToDiskPartitioned','Done writing parallel mesh',Level=8)

!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDiskPartitioned
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Generate element edge (faces in 3D) tables for given mesh.
!> Currently only for triangles and tetras. If mesh already
!> has edges do nothing.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges( Mesh, FindEdges, FindFaces )
!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     LOGICAL, OPTIONAL :: FindEdges, FindFaces

     LOGICAL :: FindEdges3D, FindFaces3d
     INTEGER :: MeshDim, SpaceDim, MaxElemDim 

     IF(PRESENT(FindEdges)) THEN
       FindEdges3D = FindEdges
     ELSE
       FindEdges3D = .TRUE.
     END IF

     IF(PRESENT(FindFaces)) THEN
       FindFaces3D = FindFaces
     ELSE
       FindFaces3D = .TRUE.
     END IF

!------------------------------------------------------------------------------

     SpaceDim = CoordinateSystemDimension()
     MeshDim = Mesh % MeshDim

     IF( MeshDim == 0 ) THEN
       CALL Fatal('FindMeshEdges','Mesh dimension is zero!')
     END IF
     IF( SpaceDim > MeshDim ) THEN
       CALL Warn('FindMeshEdges','Mesh dimension and space dimension differ: '&
           // TRIM(I2S(MeshDim))//' vs. '//TRIM(I2S(SpaceDim)))
     END IF

     MaxElemDim = EnsureElemDim( MeshDim ) 
     IF( MaxElemDim < MeshDim ) THEN
       CALL Warn('FindMeshEdges','Element dimension smaller than mesh dimension: '//&
           TRIM(I2S(MaxElemDim))//' vs '//TRIM(I2S(MeshDim)))
     END IF


     SELECT CASE( MaxElemDim )

     CASE(2)
       IF ( .NOT.ASSOCIATED( Mesh % Edges ) ) THEN
         CALL Info('FindMeshEdges','Determining edges in 2D mesh',Level=8)
         CALL FindMeshEdges2D( Mesh )
       END IF

     CASE(3)
       IF ( .NOT.ASSOCIATED(Mesh % Faces) .AND. FindFaces3D ) THEN
         CALL Info('FindMeshEdges','Determining faces in 3D mesh',Level=8)
         CALL FindMeshFaces3D( Mesh )
       END IF
       IF(FindEdges3D) THEN
         IF ( .NOT.ASSOCIATED( Mesh % Edges) ) THEN
           CALL Info('FindMeshEdges','Determining edges in 3D mesh',Level=8)
           CALL FindMeshEdges3D( Mesh )
         END IF
       END IF
     END SELECT

     CALL AssignConstraints()

CONTAINS

  ! Check that the element dimension really follows the mesh dimension
  ! The default is the MeshDim so we return immediately after that is 
  ! confirmed. 
  !--------------------------------------------------------------------
    FUNCTION EnsureElemDim(MeshDim) RESULT (MaxElemDim)

      INTEGER :: MeshDim, MaxElemDim 
      INTEGER :: i,ElemDim, ElemCode

      MaxElemDim = 0

      DO i=1,Mesh % NumberOfBulkElements
        ElemCode = Mesh % Elements(i) % Type % ElementCode
        IF( ElemCode > 500 ) THEN
          ElemDim = 3 
        ELSE IF( ElemCode > 300 ) THEN
          ElemDim = 2
        ELSE IF( ElemCode > 200 ) THEN
          ElemDim = 1
        END IF
        MaxElemDim = MAX( MaxElemDim, ElemDim ) 
        IF( MaxElemDim == MeshDim ) EXIT
      END DO
          
    END FUNCTION EnsureElemDim


    SUBROUTINE AssignConstraints()

      INTEGER, POINTER :: FaceInd(:)
      INTEGER :: i,j,k,l,n,nd,nfound
      TYPE(Element_t), POINTER :: Element, Boundary, Face, Faces(:)

      DO i=1,Mesh % NumberOfBoundaryElements
        Boundary => Mesh % Elements(Mesh % NumberOfBulkElements+i)

        Element  => Boundary % BoundaryInfo % Left
        IF (.NOT.ASSOCIATED(Element) ) &
          Element  => Boundary % BoundaryInfo % Right
        IF (.NOT.ASSOCIATED(Element) ) CYCLE

        SELECT CASE(Boundary % TYPE % DIMENSION)
        CASE(1)
          nd = Element % TYPE % NumberOfEdges
          Faces   => Mesh % Edges
          FaceInd => Element % EdgeIndexes
        CASE(2)
          nd = Element % TYPE % NumberOfFaces
          Faces   => Mesh % Faces
          FaceInd => Element % FaceIndexes
        CASE DEFAULT
          Faces => NULL()
          FaceInd => NULL()
        END SELECT

        IF ( .NOT. ASSOCIATED(Faces) .OR. .NOT. ASSOCIATED(FaceInd) ) CYCLE

        DO j=1,nd
          IF(FaceInd(j)<=0) CYCLE

          Face => Faces(FaceInd(j))
          IF ( .NOT.ASSOCIATED(Face % TYPE,Boundary % TYPE) ) CYCLE

          n = Boundary % TYPE % NumberOfNodes
          nfound = 0
          DO k=1,n
            DO l=1,n
              IF ( Boundary % NodeIndexes(k)==Face % NodeIndexes(l) ) &
                nfound = nfound+1
            END DO
          END DO
          IF ( nfound==n ) THEN
            Face % BoundaryInfo % Constraint = Boundary % BoundaryInfo % Constraint; EXIT
          END IF
        END DO
      END DO
    END SUBROUTINE AssignConstraints
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Find 2D mesh edges.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges2D( Mesh, BulkMask )
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: BulkMask(:)
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node,Edge
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
     
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    TYPE(Element_t), POINTER :: Element, Edges(:)

    LOGICAL :: Found,Masked, LG
    INTEGER :: i,j,k,n,NofEdges,Edge,Swap,Node1,Node2,istat,Degree,maxedges,allocstat
!------------------------------------------------------------------------------
!
!   Initialize:
!   -----------

    CALL Info('FindMeshEdges2D','Finding mesh edges in 2D mesh',Level=12)
    
    Masked = PRESENT(BulkMask)
    
    DO i=1,SIZE(Mesh % Elements)
       Element => Mesh % Elements(i)

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)

           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
               j=Element % Boundaryinfo % Right % ElementIndex
           END IF

           IF(j==-1) CYCLE
         END IF
         IF ( .NOT.BulkMask(j)) CYCLE
       END IF

       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector( Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    CALL Info('FindMeshEdges2D','Creating hash table of size '&
        //TRIM(I2S(Mesh % NumberOfNodes))//' for node-to-node connectivity',Level=20)
    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
      NULLIFY( HashTable(i) % Head )
    END DO
    CALL Info('FindMeshEdges2D','Hash table allocated',Level=25)
     
!------------------------------------------------------------------------------

#if 1
    Edges => NULL()
    NofEdges = 0
1   DO i=1,SIZE(Mesh % Elements)

       Element => Mesh % Elements(i)

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
               j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)
           
           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
                 j=Element % Boundaryinfo % Right % ElementIndex
           END IF
           
           IF(j==-1) CYCLE
         END IF
         
         IF(.NOT. BulkMask(j)) CYCLE
       END IF

       SELECT CASE( Element % TYPE % ElementCode / 100 )
       CASE(1) 
         CYCLE
       CASE(2)
         n = 1
       CASE(3)
         n = 3
       CASE(4)
         n = 4
       END SELECT
       
!      Loop over every edge of every element:
!      --------------------------------------
       DO k=1,n
!         We use MIN(Node1,Node2) as the hash table key:
!         ----------------------------------------------
         Node1 = Element % NodeIndexes(k)
         IF(n==1) THEN
           Node2 = Element % NodeIndexes(2)
         ELSE IF ( k<n ) THEN
           Node2 = Element % NodeIndexes(k+1)
         ELSE
           Node2 = Element % NodeIndexes(1)
         END IF
         
         IF ( Node2 < Node1 ) THEN
           Swap  = Node1
           Node1 = Node2
           Node2 = Swap
         END IF
         
!         Look the edge from the hash table:
!         ----------------------------------
         HashPtr => HashTable(Node1) % Head
         Found = .FALSE.         
         DO WHILE( ASSOCIATED( HashPtr ) )
           IF ( HashPtr % Node == Node2 ) THEN
             Found = .TRUE.
             Edge = HashPtr % Edge
             EXIT
           END IF
           HashPtr => HashPtr % Next
         END DO

         IF(.NOT. ASSOCIATED( Edges ) ) THEN
           ! Edge has already been numbered
           IF(Found ) CYCLE

           ! This is visited only the first round when Edges have not been allocated.           
           NofEdges = NofEdges + 1
           Edge = NofEdges
           
           ! Update the hash table:
           !----------------------
           ALLOCATE( HashPtr, STAT=allocstat )
           IF( allocstat /= 0 ) THEN
             CALL Fatal('FindMeshEdges2D','Allocation error for HashPtr allocation')
           END IF           
           HashPtr % Edge = Edge
           HashPtr % Node = Node2
           HashPtr % Next => HashTable(Node1) % Head
           HashTable(Node1) % Head => HashPtr
         
         ELSE 
           IF(.NOT. Found ) THEN
             CALL Fatal('FindMeshEdges2D','We should find the edge in the hash table!')
           END IF
           IF( Edge > SIZE( Edges ) ) THEN
             CALL Fatal('FindMeshEdges2D','Number of edges larger than expected!')
           END IF
                      
           IF(.NOT. ASSOCIATED(Edges(Edge) % TYPE ) ) THEN
             Degree = Element % TYPE % BasisFunctionDegree             

             Edges(Edge) % ElementIndex = Edge
             CALL AllocateVector( Edges(Edge) % NodeIndexes, Degree+1)
             ALLOCATE( Edges(Edge) % BoundaryInfo, STAT=allocstat )
             IF( allocstat /= 0 ) THEN
               CALL Fatal('FindMeshEdges2D','Allocation error for BoyndaryInfo allocation')
             END IF
             Edges(Edge) % TYPE => GetElementType( 201+Degree, .FALSE. )

             Edges(Edge) % NodeIndexes(1) = Element % NodeIndexes(k)
             IF ( k < n ) THEN
               Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(k+1)
             ELSE
               Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(1)
             END IF

             DO j=2,Degree
               Edges(Edge) % NodeIndexes(j+1) = Element % NodeIndexes(k+n+j-2)
             END DO
             
             ! Create P element definitions if needed
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
               CALL AllocatePDefinitions(Edges(Edge))
               Edges(Edge) % PDefs % P = 0
             ELSE
               NULLIFY( Edges(Edge) % PDefs )
             END IF

             Edges(Edge) % NDofs = 0
             IF (Element % NDOFs /= 0 ) Edges(Edge) % NDOFs = &
                 Element % NDOFs / Element % TYPE % NumberOfNodes * &
                 Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             NULLIFY( Edges(Edge) % EdgeIndexes )
             NULLIFY( Edges(Edge) % FaceIndexes )
             
             Edges(Edge) % BoundaryInfo % Left  => NULL()
             Edges(Edge) % BoundaryInfo % Right => NULL()
           END IF

           ! These stuctures need to be updated to both new and old edge.
           Element % EdgeIndexes(k) = Edge
           IF (i <= Mesh % NumberofBulkElements) THEN
             IF(ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
               Edges(Edge) % BoundaryInfo % Right => Element
             ELSE
               Edges(Edge) % BoundaryInfo % Left => Element
             END IF
           END IF
           
         END IF
       END DO
     END DO

     IF(.NOT. ASSOCIATED( Edges ) ) THEN
       CALL Info('FindMeshEdges2D','Allocating edge table of size: '//TRIM(I2S(NofEdges)),Level=12)
       CALL AllocateVector( Mesh % Edges, NofEdges ) 
       Edges => Mesh % Edges
       GOTO 1
     END IF
     
#else
     maxedges = 0
     DO i=1,SIZE( Mesh % Elements )
       Element => Mesh % Elements(i)
       maxedges = MAX( maxedges, Element % TYPE % NumberOfEdges )
     END DO
              
     n = maxedges * Mesh % NumberOfBulkElements
     CALL Info('FindMeshEdges2D','Allocating edge table of size: '//TRIM(I2S(n)),Level=12)
     CALL AllocateVector( Mesh % Edges, n ) 
     Edges => Mesh % Edges
     
!   Loop over elements:
!   -------------------
    NofEdges = 0
    DO i=1,SIZE(Mesh % Elements)

       Element => Mesh % Elements(i)

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)

           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
               j=Element % Boundaryinfo % Right % ElementIndex
           END IF

           IF(j==-1) CYCLE
         END IF

         IF(.NOT. BulkMask(j)) CYCLE
       END IF


       SELECT CASE( Element % TYPE % ElementCode / 100 )
         CASE(1) 
            CYCLE
         CASE(2)
            n = 1
         CASE(3)
            n = 3
         CASE(4)
            n = 4
       END SELECT

!      Loop over every edge of every element:
!      --------------------------------------
       DO k=1,n
!         We use MIN(Node1,Node2) as the hash table key:
!         ----------------------------------------------
          Node1 = Element % NodeIndexes(k)
          IF(n==1) THEN
             Node2 = Element % NodeIndexes(2)
          ELSE IF ( k<n ) THEN
             Node2 = Element % NodeIndexes(k+1)
          ELSE
             Node2 = Element % NodeIndexes(1)
          END IF

          IF ( Node2 < Node1 ) THEN
             Swap  = Node1
             Node1 = Node2
             Node2 = Swap
          END IF

!         Look the edge from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.         
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node == Node2 ) THEN
                Found = .TRUE.
                Edge = HashPtr % Edge
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO

!         Existing edge, update structures:
!         ----------------------------------
          IF ( Found ) THEN
             Element % EdgeIndexes(k) = Edge
             IF (i<=Mesh % NumberofBulkElements) THEN
               IF(ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
                 Edges(Edge) % BoundaryInfo % Right => Element
               ELSE
                 Edges(Edge) % BoundaryInfo % Left => Element
               END IF
             END IF
          ELSE

!            Edge not yet there, create:
!            ---------------------------
             NofEdges = NofEdges + 1
             Edge = NofEdges

             Degree = Element % TYPE % BasisFunctionDegree
             
             IF( Edge > SIZE( Edges ) ) THEN
               CALL Fatal('FindMeshEdges2D','Number of edges larger than expected!')
             END IF

             Edges(Edge) % ElementIndex = Edge
             CALL AllocateVector( Edges(Edge) % NodeIndexes, Degree+1)
             ALLOCATE( Edges(Edge) % BoundaryInfo, STAT=allocstat )
             IF( allocstat /= 0 ) THEN
               CALL Fatal('FindMeshEdges2D','Allocation error for BoyndaryInfo allocation')
             END IF

             Edges(Edge) % TYPE => GetElementType( 201+Degree, .FALSE. )

             Edges(Edge) % NodeIndexes(1) = Element % NodeIndexes(k)
             IF ( k < n ) THEN
                Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(k+1)
             ELSE
                Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(1)
             END IF

             DO j=2,Degree
                Edges(Edge) % NodeIndexes(j+1) = Element % NodeIndexes(k+n+j-2)
             END DO
             
             ! Create P element definitions if needed
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
               CALL AllocatePDefinitions(Edges(Edge))
               Edges(Edge) % PDefs % P = 0
             ELSE
               NULLIFY( Edges(Edge) % PDefs )
             END IF

             Edges(Edge) % NDofs = 0
             IF (Element % NDOFs /= 0 ) Edges(Edge) % NDOFs = &
                 Element % NDOFs / Element % TYPE % NumberOfNodes * &
                 Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             NULLIFY( Edges(Edge) % EdgeIndexes )
             NULLIFY( Edges(Edge) % FaceIndexes )

             Element % EdgeIndexes(k) = Edge
             Edges(Edge) % BoundaryInfo % Left  => Null()
             Edges(Edge) % BoundaryInfo % Right => Null()
             IF(i<=Mesh % NumberOfBulkElements) Edges(Edge) % BoundaryInfo % Left => Element
              
!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr, STAT=allocstat )
             IF( allocstat /= 0 ) THEN
               CALL Fatal('FindMeshEdges2D','Allocation error for HashPtr allocation')
             END IF

             HashPtr % Edge = Edge
             HashPtr % Node = Node2
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO
#endif
    
    Mesh % NumberOfEdges = NofEdges
    CALL Info('FindMeshEdges2D','Number of edges found: '//TRIM(I2S(NofEdges)),Level=10)

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )

    CALL Info('FindMeshEdges2D','All done',Level=12)

!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh faces.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshFaces3D( Mesh, BulkMask)
    USE PElementMaps, ONLY : GetElementFaceMap
    USE PElementBase, ONLY : isPTetra

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: BulkMask(:)
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node1,Node2,Face
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
    
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    LOGICAL :: Found,Masked,LG
    INTEGER :: n1,n2,n3,n4
    INTEGER :: i,j,k,n,NofFaces,Face,Swap,Node1,Node2,Node3,istat,Degree,facenodes
     
    TYPE(Element_t), POINTER :: Element, Faces(:)

    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET  :: TetraFaceMap(4,6), BrickFaceMap(6,9), &
         WedgeFaceMap(5,8), PyramidFaceMap(5,8), TriFaceMap(1,3), QuadFaceMap(1,4)
    
    INTEGER :: nf(4)
!------------------------------------------------------------------------------
    
    CALL Info('FindMeshFaces3D','Finding mesh faces in 3D mesh',Level=12)

    Masked = PRESENT(BulkMask)

    TriFaceMap(1,:)  = [1,2,3]
    QuadFaceMap(1,:) = [1,2,3,4]

    TetraFaceMap(1,:) = [ 1, 2, 3, 5, 6, 7 ]
    TetraFaceMap(2,:) = [ 1, 2, 4, 5, 9, 8 ]
    TetraFaceMap(3,:) = [ 2, 3, 4, 6, 10, 9 ]
    TetraFaceMap(4,:) = [ 3, 1, 4, 7, 8,10 ]

    WedgeFaceMap(1,:) = [ 1, 2, 3, 7, 8, 9, -1, -1 ]
    WedgeFaceMap(2,:) = [ 4, 5, 6, 10, 11, 12, -1, -1 ]
    WedgeFaceMap(3,:) = [ 1, 2, 5, 4, 7, 14, 10, 13 ]
    WedgeFaceMap(4,:) = [ 3, 2, 5, 6, 8, 14, 11, 15 ]
    WedgeFaceMap(5,:) = [ 3, 1, 4, 6, 9, 13, 12, 15 ]

    PyramidFaceMap(1,:) = [ 1, 2, 3, 4,  6,  7,  8,  9 ]
    PyramidFaceMap(2,:) = [ 1, 2, 5, 6, 11, 10, -1, -1 ]
    PyramidFaceMap(3,:) = [ 2, 3, 5, 7, 12, 11, -1, -1 ]
    PyramidFaceMap(4,:) = [ 3, 4, 5, 8, 13, 12, -1, -1 ]
    PyramidFaceMap(5,:) = [ 4, 1, 5, 9, 10, 13, -1, -1 ]

    BrickFaceMap(1,:) = [ 1, 2, 3, 4,  9, 10, 11, 12, 25 ]
    BrickFaceMap(2,:) = [ 5, 6, 7, 8, 17, 18, 19, 20, 26 ]
    BrickFaceMap(3,:) = [ 1, 2, 6, 5,  9, 14, 17, 13, 21 ]
    BrickFaceMap(4,:) = [ 2, 3, 7, 6, 10, 15, 18, 14, 22 ]
    BrickFaceMap(5,:) = [ 3, 4, 8, 7, 11, 16, 19, 15, 23 ]
    BrickFaceMap(6,:) = [ 4, 1, 5, 8, 12, 13, 20, 16, 24 ]

!
!   Initialize:
!   -----------   
    DO i=1,SIZE(Mesh % Elements)
       Element => Mesh % Elements(i)

       IF(.NOT.ASSOCIATED(Element % Type)) CYCLE
       IF(Element % Type % ElementCode<300 ) CYCLE

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)

           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
               j=Element % Boundaryinfo % Right % ElementIndex
           END IF

           IF(j==-1) CYCLE
         END IF

         IF(.NOT. BulkMask(j)) CYCLE
       END IF

       IF ( .NOT. ASSOCIATED( Element % FaceIndexes ) ) &
          CALL AllocateVector(Element % FaceIndexes, Element % TYPE % NumberOfFaces )
       Element % FaceIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

#if 1
!   Loop over elements:
!   -------------------
    NofFaces = 0
    Faces => NULL()

1   DO i=1,SIZE(Mesh % Elements)

      Element => Mesh % Elements(i)
      IF(.NOT.ASSOCIATED(Element % Type)) CYCLE
      IF(Element % Type % ElementCode < 300 ) Cycle

      IF(Masked) THEN
        j = i
        IF(i>Mesh % NumberOfBulkElements) THEN
          j = -1
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

          LG=.FALSE.
          IF(j>0) LG=BulkMask(j)

          IF(.NOT. LG) THEN
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
                j=Element % Boundaryinfo % Right % ElementIndex
          END IF

          IF(j==-1) CYCLE
        END IF
        IF(.NOT. BulkMask(j)) CYCLE
      END IF

      ! For P elements mappings are different
      IF ( ASSOCIATED(Element % PDefs) ) THEN
        CALL GetElementFaceMap(Element, FaceMap)
        n = Element % TYPE % NumberOfFaces
      ELSE
        SELECT CASE( Element % TYPE % ElementCode / 100 )
        CASE(3)
          n = 1
          FaceMap => TriFaceMap
        CASE(4)
          n = 1
          FaceMap => QuadFaceMap
        CASE(5)
          n = 4
          FaceMap => TetraFaceMap
        CASE(6)
          n = 5
          FaceMap => PyramidFaceMap
        CASE(7)
          n = 5 
          FaceMap => WedgeFaceMap
        CASE(8)
          n = 6
          FaceMap => BrickFaceMap
        CASE DEFAULT
          CALL Fatal('FindMeshFaces','Element type '&
              //TRIM(I2S(Element % Type % ElementCode))//' not implemented!')
        END SELECT
      END IF
 
!      Loop over every face of every element:
!      --------------------------------------
      DO k=1,n
                    
        SELECT CASE( Element % TYPE % ElementCode / 100 )
          
        CASE(3)
          ! Triangle:
          !=======
          facenodes = 3

        CASE(4)
          ! Quad:
          !=======
          facenodes = 4

        CASE(5)
          ! Tetras:
          !=======
          facenodes = 3

        CASE(6)
          ! Pyramids:
          !=========
          IF ( k == 1 ) THEN
            facenodes = 4
          ELSE
            facenodes = 3
          END IF
          
        CASE(7)
          ! Wedges:
          !=======
          IF ( k <= 2 ) THEN
            facenodes = 3
          ELSE
            facenodes = 4
          END IF
                
        CASE(8)
          ! Bricks:
          !=======
          facenodes = 4
          
        CASE DEFAULT
          WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
          CALL Fatal('FindMeshFaces',Message)
        END SELECT

        nf(1:facenodes) = Element % NodeIndexes(FaceMap(k,1:facenodes))
        CALL sort( facenodes, nf )
        
!         We use MIN(Node1,Node2,Node3) as the hash table key:
!         ---------------------------------------------------
        Node1 = nf(1)
        Node2 = nf(2)
        Node3 = nf(3)
          
!         Look the face from the hash table:
!         ----------------------------------
        HashPtr => HashTable(Node1) % Head
        Found = .FALSE.
        DO WHILE( ASSOCIATED( HashPtr ) )
          IF ( HashPtr % Node1 == Node2 .AND. HashPtr % Node2 == Node3) THEN
            Found = .TRUE.
            Face = HashPtr % Face
            EXIT
          END IF
          HashPtr => HashPtr % Next
        END DO
        
!         Existing face, update structures:
!         ----------------------------------

        IF( .NOT. ASSOCIATED( Faces ) ) THEN
          IF(Found ) CYCLE

          ! Update the hash table:
          !----------------------
          NofFaces = NofFaces + 1
          Face = NofFaces
          ALLOCATE( HashPtr )
          HashPtr % Face = Face
          HashPtr % Node1 = Node2
          HashPtr % Node2 = Node3
          HashPtr % Next => HashTable(Node1) % Head
          HashTable(Node1) % Head => HashPtr
        ELSE
          IF(.NOT. Found ) THEN
            CALL Fatal('FindMeshFaces3D','We should find the edge in the hash table!')
          END IF
          IF( Face > SIZE( Faces ) ) THEN
            CALL Fatal('FindMeshFaces3D','Number of faces larger than expected!')
          END IF
          
          IF(.NOT. ASSOCIATED( Faces(Face) % TYPE ) ) THEN
            ! Face not yet there, create:
            !---------------------------
            Degree = Element % TYPE % BasisFunctionDegree
            Faces(Face) % ElementIndex = Face
            
            SELECT CASE( Element % TYPE % ElementCode / 100 )

            CASE(1,2)
              CYCLE

            CASE(3)
              ! linear tri
              !-----------
              SELECT CASE( Degree ) 
              CASE(1)
                n1 = 3
              CASE DEFAULT
              END SELECT
              
              Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              
            CASE(4)
              ! linear quad
              !-----------
              SELECT CASE( Degree ) 
              CASE(1)
                n1 = 4
              CASE DEFAULT
              END SELECT              
              Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
              
            CASE(5)
              ! for tetras:
              !-----------
              SELECT CASE( Degree ) 
              CASE(1)
                n1 = 3
              CASE(2)
                n1 = 6
              CASE(3)
                n1 = 10
              END SELECT
              
              Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              
            CASE(6)              
               ! Pyramids ( 605 and 613 supported )
               !-------------------------------
              IF ( k == 1 ) THEN
                n1 = Degree * 4
                Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
              ELSE
                n1 = Degree * 3
                Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              END IF
              
            CASE(7)
               ! for wedges, 706 and 715 supported:
               !-------------------------------
              IF ( k <= 2 ) THEN
                n1 = Degree * 3
                Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              ELSE
                n1 = Degree * 4
                Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
              END IF
              
            CASE(8)
               ! for bricks:
               !-----------
              SELECT CASE( Element % TYPE % NumberOfNodes ) 
              CASE(8)
                n1 = 4
              CASE(20)
                n1 = 8
              CASE(27)
                n1 = 9
              END SELECT
              
              Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE.)
              
            CASE DEFAULT
              CALL Fatal('FindMeshFaces','Element type '&
                  //TRIM(I2S(Element % TYPE % ElementCode))//' not implemented!')
              
            END SELECT
            
             ! Allocate p structures for p elements
            IF ( ASSOCIATED( Element % PDefs ) ) THEN
              CALL AllocatePDefinitions(Faces(Face))
              Faces(Face) % PDefs % P = 0
            ELSE
              NULLIFY( Faces(Face) % PDefs )
            END IF
            
            Faces(Face) % NDOFs  = 0
            IF (Element % NDOFs /= 0) Faces(Face) % NDOFs = &
                Element % NDOFs / Element % TYPE % NumberOfNodes * &
                Faces(Face) % TYPE % NumberOfNodes
            Faces(Face) % BDOFs  = 0
            Faces(Face) % DGDOFs = 0
            Faces(Face) % EdgeIndexes => NULL()
            Faces(Face) % FaceIndexes => NULL()
            
            CALL AllocateVector( Faces(Face) % NodeIndexes,n1 )
            DO n2=1,n1
              Faces(Face) % NodeIndexes(n2) = &
                  Element % NodeIndexes(FaceMap(k,n2)) 
            END DO
            
            ALLOCATE( Faces(Face) % BoundaryInfo )
            Faces(Face) % BoundaryInfo % Left  => NULL()
            Faces(Face) % BoundaryInfo % Right => NULL()
          END IF

          Element % FaceIndexes(k) = Face            
          IF(i<=Mesh % NumberOfBulkElements) THEN
            IF( ASSOCIATED(Faces(Face) % BoundaryInfo % Left) ) THEN
              Faces(Face) % BoundaryInfo % Right => Element
            ELSE
              Faces(Face) % BoundaryInfo % Left => Element
            END IF
          END IF
          
        END IF
      END DO
    END DO

    IF(.NOT. ASSOCIATED( Faces ) ) THEN
      CALL Info('FindMeshFaces3D','Allocating face table of size: '&
          //TRIM(I2S(NofFaces)),Level=25)
      CALL AllocateVector( Mesh % Faces, NofFaces, 'FindMeshFaces3D' )
      Faces => Mesh % Faces
      GOTO 1
    END IF
    
#else
   
    IF(Masked) THEN
      n = 6 * COUNT(BulkMask)
    ELSE
      n = 6 * Mesh % NumberOfBulkElements
    END IF

    CALL Info('FindMeshFaces3D','Allocating face table of size: '&
        //TRIM(I2S(n)),Level=25)
    CALL AllocateVector( Mesh % Faces, n, 'FindMeshFaces3D' )
    Faces => Mesh % Faces

    
!   Loop over elements:
!   -------------------
    NofFaces = 0
    DO i=1,SIZE(Mesh % Elements)
 
       Element => Mesh % Elements(i)
       IF(.NOT.ASSOCIATED(Element % Type)) CYCLE
       IF(Element % Type % ElementCode<300 ) Cycle

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)

           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
               j=Element % Boundaryinfo % Right % ElementIndex
           END IF

           IF(j==-1) CYCLE
         END IF
         IF(.NOT. BulkMask(j)) CYCLE
       END IF


       ! For P elements mappings are different
       IF ( ASSOCIATED(Element % PDefs) ) THEN
          CALL GetElementFaceMap(Element, FaceMap)
          n = Element % TYPE % NumberOfFaces
       ELSE
          SELECT CASE( Element % TYPE % ElementCode / 100 )
          CASE(3)
             n = 1
             FaceMap => TriFaceMap
          CASE(4)
             n = 1
             FaceMap => QuadFaceMap
          CASE(5)
             n = 4
             FaceMap => TetraFaceMap
          CASE(6)
             n = 5
             FaceMap => PyramidFaceMap
          CASE(7)
             n = 5 
             FaceMap => WedgeFaceMap
          CASE(8)
             n = 6
             FaceMap => BrickFaceMap
          CASE DEFAULT
             CYCLE
             ! WRITE(Message,*) 'Element type',Element % Type % ElementCode,'not implemented.' 
             ! CALL Fatal('FindMeshFaces',Message)
          END SELECT
       END IF
 
!      Loop over every face of every element:
!      --------------------------------------
       DO k=1,n
          
          
!         We use MIN(Node1,Node2,Node3) as the hash table key:
!         ---------------------------------------------------
          SELECT CASE( Element % TYPE % ElementCode / 100 )
             CASE(3)
!
!               Ttriangle:
!               =======
                nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                CALL sort( 3, nf )
             CASE(4)
!
!               Quad:
!               =======
                nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                CALL sort( 4, nf )

             CASE(5)
!
!               Tetras:
!               =======
                nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                CALL sort( 3, nf )

             CASE(6)
!
!               Pyramids:
!               =========
                IF ( k == 1 ) THEN
                   nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                   CALL sort( 4, nf )
                ELSE
                   nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                   CALL sort( 3, nf )
                END IF

             CASE(7)
!
!               Wedges:
!               =======
                IF ( k <= 2 ) THEN
                   nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                   CALL sort( 3, nf )
                ELSE
                   nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                   CALL sort( 4, nf )
                END IF
                
             CASE(8)
!
!               Bricks:
!               =======
                nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                CALL sort( 4, nf )

             CASE DEFAULT
                WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
                CALL Fatal('FindMeshFaces',Message)
          END SELECT

          Node1 = nf(1)
          Node2 = nf(2)
          Node3 = nf(3)
          
!         Look the face from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node1 == Node2 .AND. HashPtr % Node2 == Node3) THEN
                Found = .TRUE.
                Face = HashPtr % Face
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO

!         Existing face, update structures:
!         ----------------------------------
          IF ( Found ) THEN
             Element % FaceIndexes(k) = Face
             IF(i<=Mesh % NumberOfBulkElements) THEN
               IF( ASSOCIATED(Faces(Face) % BoundaryInfo % Left) ) THEN
                 Faces(Face) % BoundaryInfo % Right => Element
               ELSE
                 Faces(Face) % BoundaryInfo % Left => Element
               END IF
             END IF
          ELSE

!            Face not yet there, create:
!            ---------------------------
             NofFaces = NofFaces + 1
             Face = NofFaces
             Faces(Face) % ElementIndex = Face

             Degree = Element % TYPE % BasisFunctionDegree


             SELECT CASE( Element % TYPE % ElementCode / 100 )
             CASE(3)
               !
               !               linear tri
               !               -----------
               SELECT CASE( Degree ) 
               CASE(1)
                 n1 = 3
               CASE DEFAULT
               END SELECT

               Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )

             CASE(4)
               !
               !               linear quad
               !               -----------
               SELECT CASE( Degree ) 
               CASE(1)
                 n1 = 4
               CASE DEFAULT
               END SELECT

               Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )

             CASE(5)
               !
               !               for tetras:
               !               -----------
               SELECT CASE( Degree ) 
               CASE(1)
                 n1 = 3
               CASE(2)
                 n1 = 6
               CASE(3)
                 n1 = 10
               END SELECT

               Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )

             CASE(6)

               !               Pyramids ( 605 and 613 supported )
               !               -------------------------------
               IF ( k == 1 ) THEN
                 n1 = Degree * 4
                 Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
               ELSE
                 n1 = Degree * 3
                 Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
               END IF

             CASE(7)

               !               for wedges, 706 and 715 supported:
               !               -------------------------------
               IF ( k <= 2 ) THEN
                 n1 = Degree * 3
                 Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
               ELSE
                 n1 = Degree * 4
                 Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
               END IF


             CASE(8)
               !
               !               for bricks:
               !               -----------
               SELECT CASE( Element % TYPE % NumberOfNodes ) 
               CASE(8)
                 n1 = 4
               CASE(20)
                 n1 = 8
               CASE(27)
                 n1 = 9
               END SELECT

               Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE.)

             CASE DEFAULT
               WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
               CALL Fatal('FindMeshFaces',Message)

             END SELECT

             ! Allocate p structures for p elements
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
                CALL AllocatePDefinitions(Faces(Face))
                Faces(Face) % PDefs % P = 0
             ELSE
               NULLIFY( Faces(Face) % PDefs )
             END IF
             
             Faces(Face) % NDOFs  = 0
             IF (Element % NDOFs /= 0) Faces(Face) % NDOFs = &
                 Element % NDOFs / Element % TYPE % NumberOfNodes * &
                      Faces(Face) % Type % NumberOfNodes
             Faces(Face) % BDOFs  = 0
             Faces(Face) % DGDOFs = 0
             Faces(Face) % EdgeIndexes => NULL()
             Faces(Face) % FaceIndexes => NULL()

             CALL AllocateVector( Faces(Face) % NodeIndexes,n1 )
             DO n2=1,n1
                Faces(Face) % NodeIndexes(n2) = &
                         Element % NodeIndexes(FaceMap(k,n2)) 
             END DO

             Element % FaceIndexes(k) = Face

             ALLOCATE( Faces(Face) % BoundaryInfo )
             Faces(Face) % BoundaryInfo % Left  => Null()
             Faces(Face) % BoundaryInfo % Right => Null()
             IF(i<=Mesh % NumberOfBulkElements) Faces(Face) % BoundaryInfo % Left => Element
              
!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr )
             HashPtr % Face = Face
             HashPtr % Node1 = Node2
             HashPtr % Node2 = Node3
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO
#endif
    
    Mesh % NumberOfFaces = NofFaces
    CALL Info('FindMeshFaces3D','Number of faces found: '//TRIM(I2S(NofFaces)),Level=10)

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )

    CALL Info('FindMeshFaces3D','All done',Level=12)
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshFaces3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh edges.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges3D( Mesh )
    USE PElementMaps, ONLY : GetElementEdgeMap, GetElementFaceEdgeMap
    USE PElementBase, ONLY : isPPyramid

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node1,Edge
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
    
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    LOGICAL :: Found
    INTEGER :: n1,n2, n_e, maxedges
    INTEGER :: i,j,k,n,NofEdges,Edge,Node1,Node2,istat,Degree,ii,jj
     
    TYPE(Element_t), POINTER :: Element, Edges(:), Face

    INTEGER, POINTER :: EdgeMap(:,:), FaceEdgeMap(:,:)
    INTEGER, TARGET  :: TetraEdgeMap(6,3), BrickEdgeMap(12,3), TetraFaceMap(4,6), &
      WedgeEdgeMap(9,3), PyramidEdgeMap(8,3), TetraFaceEdgeMap(4,3), &
      BrickFaceEdgeMap(8,4), WedgeFaceEdgeMap(6,4), PyramidFaceEdgeMap(5,4), &
         QuadEdgeMap(4,3), TriEdgeMap(3,3), TriFaceMap(1,3), QuadFaceMap(1,4), LineEdgeMap(1,2)
!------------------------------------------------------------------------------
    
    CALL Info('FindMeshEdges3D','Finding mesh edges in 3D mesh',Level=12)

    LineEdgeMap(1,:) = [1,2]

    TriEdgeMap(1,:) = [1,2,4]
    TriEdgeMap(2,:) = [2,3,5]
    TriEdgeMap(3,:) = [3,1,6]

    TriFaceMap(1,:) = [1,2,3]

    QuadEdgeMap(1,:) = [1,2,5]
    QuadEdgeMap(2,:) = [2,3,6]
    QuadEdgeMap(3,:) = [3,4,7]
    QuadEdgeMap(4,:) = [4,1,8]

    QuadFaceMap(1,:) = [1,2,3,4]

    TetraFaceMap(1,:) = [ 1, 2, 3, 5, 6, 7 ]
    TetraFaceMap(2,:) = [ 1, 2, 4, 5, 9, 8 ]
    TetraFaceMap(3,:) = [ 2, 3, 4, 6,10, 9 ]
    TetraFaceMap(4,:) = [ 3, 1, 4, 7, 8,10 ]

    TetraFaceEdgeMap(1,:) = [ 1,2,3 ]
    TetraFaceEdgeMap(2,:) = [ 1,5,4 ]
    TetraFaceEdgeMap(3,:) = [ 2,6,5 ]
    TetraFaceEdgeMap(4,:) = [ 3,4,6 ]

    TetraEdgeMap(1,:) = [ 1,2,5 ]
    TetraEdgeMap(2,:) = [ 2,3,6 ]
    TetraEdgeMap(3,:) = [ 3,1,7 ]
    TetraEdgeMap(4,:) = [ 1,4,8 ]
    TetraEdgeMap(5,:) = [ 2,4,9 ]
    TetraEdgeMap(6,:) = [ 3,4,10 ]

    PyramidEdgeMap(1,:) = [ 1,2,1 ]
    PyramidEdgeMap(2,:) = [ 2,3,1 ]
    PyramidEdgeMap(3,:) = [ 3,4,1 ]
    PyramidEdgeMap(4,:) = [ 4,1,1 ]
    PyramidEdgeMap(5,:) = [ 1,5,1 ]
    PyramidEdgeMap(6,:) = [ 2,5,1 ]
    PyramidEdgeMap(7,:) = [ 3,5,1 ]
    PyramidEdgeMap(8,:) = [ 4,5,1 ]

    PyramidFaceEdgeMap(1,:) = [ 1,2,3,4 ]
    PyramidFaceEdgeMap(2,:) = [ 1,6,5,0 ]
    PyramidFaceEdgeMap(3,:) = [ 2,7,6,0 ]
    PyramidFaceEdgeMap(4,:) = [ 3,8,7,0 ]
    PyramidFaceEdgeMap(5,:) = [ 4,5,8,0 ]

    WedgeEdgeMap(1,:) = [ 1, 2, 1 ]
    WedgeEdgeMap(2,:) = [ 2, 3, 1 ]
    WedgeEdgeMap(3,:) = [ 1, 3, 1 ]
    WedgeEdgeMap(4,:) = [ 4, 5, 1 ]
    WedgeEdgeMap(5,:) = [ 5, 6, 1 ]
    WedgeEdgeMap(6,:) = [ 6, 4, 1 ]
    WedgeEdgeMap(7,:) = [ 1, 4, 1 ]
    WedgeEdgeMap(8,:) = [ 2, 5, 1 ]
    WedgeEdgeMap(9,:) = [ 3, 6, 1 ]

    WedgeFaceEdgeMap(1,:) = [ 1,2,3,0 ]
    WedgeFaceEdgeMap(2,:) = [ 4,5,6,0 ]
    WedgeFaceEdgeMap(3,:) = [ 1,8,4,7 ]
    WedgeFaceEdgeMap(4,:) = [ 2,9,5,8 ]
    WedgeFaceEdgeMap(5,:) = [ 3,7,6,9 ]

    BrickEdgeMap(1,:) = [ 1, 2,  9 ]
    BrickEdgeMap(2,:) = [ 2, 3,  10 ]
    BrickEdgeMap(3,:) = [ 4, 3,  11 ]
    BrickEdgeMap(4,:) = [ 1, 4,  12 ]
    BrickEdgeMap(5,:) = [ 5, 6,  13 ]
    BrickEdgeMap(6,:) = [ 6, 7,  14 ]
    BrickEdgeMap(7,:) = [ 8, 7,  15 ]
    BrickEdgeMap(8,:) = [ 5, 8,  16 ]
    BrickEdgeMap(9,:) = [ 1, 5,  17 ]
    BrickEdgeMap(10,:) = [ 2, 6, 18 ]
    BrickEdgeMap(11,:) = [ 3, 7, 19 ]
    BrickEdgeMap(12,:) = [ 4, 8, 20 ]

    BrickFaceEdgeMap(1,:) = [ 1,2,3,4   ]
    BrickFaceEdgeMap(2,:) = [ 5,6,7,8   ]    
    BrickFaceEdgeMap(3,:) = [ 1,10,5,9  ]
    BrickFaceEdgeMap(4,:) = [ 2,11,6,10 ]
    BrickFaceEdgeMap(5,:) = [ 3,12,7,11 ]
    BrickFaceEdgeMap(6,:) = [ 4,9,8,12  ]

!
!   Initialize:
    !   -----------
    n_e = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

    DO i=1,n_e
       Element => Mesh % Elements(i)
       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector(Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    CALL Info('FindMeshEdges3D','Hash table allocated',Level=25)

    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

#if 1
    !   Loop over elements:
    !   -------------------
    NofEdges = 0
    Edges => NULL()
    
1   DO i=1,n_e
      Element => Mesh % Elements(i)
      
      ! For P elements mappings are different
      IF ( ASSOCIATED(Element % PDefs) ) THEN
        CALL GetElementEdgeMap( Element, EdgeMap )
        CALL GetElementFaceEdgeMap( Element, FaceEdgeMap ) 
        n = Element % TYPE % NumberOfEdges
      ELSE 
        SELECT CASE( Element % TYPE % ElementCode / 100 )
        CASE(1)
          CYCLE
        CASE(2)
          n = 1
          EdgeMap => LineEdgeMap
          FaceEdgeMap => NULL()
        CASE(3)
          n = 3
          EdgeMap => TriEdgeMap
          FaceEdgeMap => NULL()
        CASE(4)
          n = 4
          EdgeMap => QuadEdgeMap
          FaceEdgeMap => NULL()
        CASE(5)
          n = 6
          EdgeMap => TetraEdgeMap
          FaceEdgeMap => TetraFaceEdgeMap
        CASE(6)
          n = 8
          EdgeMap => PyramidEdgeMap
          FaceEdgeMap => PyramidFaceEdgeMap
        CASE(7)
          n = 9
          EdgeMap => WedgeEdgeMap
          FaceEdgeMap => WedgeFaceEdgeMap
        CASE(8)
          n = 12
          EdgeMap => BrickEdgeMap
          FaceEdgeMap => BrickFaceEdgeMap
        CASE DEFAULT
          CALL Fatal('FindMeshEdges3D','Element type '//TRIM(I2S(Element % TYPE % ElementCode))//' not implemented!') 
        END SELECT
      END IF

!      Loop over every edge of every element:
!      --------------------------------------
      DO k=1,n

!         Use MIN(Node1,Node2) as key to hash table:
!         ------------------------------------------
        n1 = Element % NodeIndexes(EdgeMap(k,1))
        n2 = Element % NodeIndexes(EdgeMap(k,2))
        IF ( n1 < n2 ) THEN
          Node1 = n1
          Node2 = n2
        ELSE
          Node1 = n2
          Node2 = n1
        END IF

        ! Look the edge from the hash table:
        !----------------------------------
        HashPtr => HashTable(Node1) % Head
        Found = .FALSE.
        DO WHILE( ASSOCIATED( HashPtr ) )
          IF ( HashPtr % Node1 == Node2 ) THEN
            Found = .TRUE.
            Edge = HashPtr % Edge
            EXIT
          END IF
          HashPtr => HashPtr % Next
        END DO
        
        IF(.NOT. ASSOCIATED( Edges ) ) THEN
          IF( Found ) CYCLE

          NofEdges = NofEdges + 1
          Edge = NofEdges
          
          ! Update the hash table:
          !----------------------
          ALLOCATE( HashPtr )
          HashPtr % Edge = Edge
          HashPtr % Node1 = Node2
          HashPtr % Next => HashTable(Node1) % Head
          HashTable(Node1) % Head => HashPtr
        ELSE
          IF(.NOT. Found ) THEN
            CALL Fatal('FindMeshEdges3D','We should find the edge in the hash table!')
          END IF
          IF( Edge > SIZE( Edges ) ) THEN
            CALL Fatal('FindMeshEdges3D','Number of edges larger than expected!')
          END IF
                    
          IF( ASSOCIATED( Edges(Edge) % TYPE ) ) THEN
            IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
              Edges(Edge) % BoundaryInfo % Left  => Element
            ELSE
              Edges(Edge) % BoundaryInfo % Right => Element
            END IF
          ELSE
            Degree = Element % TYPE % BasisFunctionDegree

            ! Edge is always a line segment with deg+1 nodes:
            !-----------------------------------------------
            Edges(Edge) % TYPE => GetElementType( 201 + degree, .FALSE.)

            Edges(Edge) % NDOFs  = 0
            IF (Element % NDOFs /= 0) Edges(Edge) % NDOFs = &
                Element % NDOFs / Element % TYPE % NumberOfNodes * &
                Edges(Edge) % TYPE % NumberOfNodes
            Edges(Edge) % BDOFs  = 0
            Edges(Edge) % DGDOFs = 0
            Edges(Edge) % EdgeIndexes => NULL()
            Edges(Edge) % FaceIndexes => NULL()
            
            CALL AllocateVector( Edges(Edge) % NodeIndexes, degree + 1 )
            DO n2=1,degree+1
              Edges(Edge) % NodeIndexes(n2) = &
                  Element % NodeIndexes(EdgeMap(k,n2))
            END DO
            
            ALLOCATE( Edges(Edge) % BoundaryInfo )
            Edges(Edge) % BoundaryInfo % Left  => NULL()
            Edges(Edge) % BoundaryInfo % Right => NULL()
            
            ! Allocate P element definitions 
            IF ( ASSOCIATED( Element % PDefs ) ) THEN
              CALL AllocatePDefinitions(Edges(Edge))              
              Edges(Edge) % PDefs % P = 0
              Edges(Edge) % PDefs % pyramidQuadEdge = .FALSE.
            ELSE
              NULLIFY( Edges(Edge) % PDefs )
            END IF            
          END IF

          ! Stuff for both existing and new edge
          !--------------------------------------
          Element % EdgeIndexes(k) = Edge
          
          ! Mark edge as an edge of pydamid square face 
          IF (isPPyramid(Element) .AND. k < 5) THEN
            Edges(Edge) % PDefs % pyramidQuadEdge = .TRUE.
          END IF
          
          IF ( ASSOCIATED(Mesh % Faces) .AND. ASSOCIATED(FaceEdgeMap) ) THEN
            DO ii=1,Element % TYPE % NumberOfFaces
              Face => Mesh % Faces(Element % FaceIndexes(ii))
              IF ( .NOT. ASSOCIATED(Face % EdgeIndexes) ) THEN
                ALLOCATE(Face % EdgeIndexes(Face % TYPE % NumberOfEdges))
                Face % EdgeIndexes = 0
              END IF
              DO jj=1,Face % TYPE % NumberOfEdges
                IF (FaceEdgeMap(ii,jj) == k) THEN
                  Face % EdgeIndexes(jj) = Edge
                  IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
                    Edges(Edge) % BoundaryInfo % Left => Face
                  ELSE
                    Edges(Edge) % BoundaryInfo % Right => Face
                  END IF
                  EXIT
                END IF
              END DO
            END DO
          END IF
        END IF
          
      END DO
    END DO

    IF(.NOT. ASSOCIATED( Edges ) ) THEN  
      CALL Info('FindMeshEdges3D','Allocating edge table of size: '//TRIM(I2S(NofEdges)),Level=20)
      CALL AllocateVector( Mesh % Edges, NofEdges ) 
      Edges => Mesh % Edges
      CALL Info('FindMeshEdges3D','Edge table allocated',Level=25)
      GOTO 1
    END IF
          
#else       
    maxedges = 0
    DO i=1,Mesh % NumberOfBulkElements 
      Element => Mesh % Elements(i)
      maxedges = MAX( maxedges, Element % Type % NumberOfEdges ) 
    END DO            
    n = maxedges*Mesh % NumberOfBulkElements  

    CALL Info('FindMeshEdges3D','Allocating edge table of size: '//TRIM(I2S(n)),Level=12)
    CALL AllocateVector( Mesh % Edges, n ) 
    Edges => Mesh % Edges

    CALL Info('FindMeshEdges3D','Edge table allocated',Level=25)
    
!   Loop over elements:
!   -------------------
    NofEdges = 0
    DO i=1,n_e
       Element => Mesh % Elements(i)

       ! For P elements mappings are different
       IF ( ASSOCIATED(Element % PDefs) ) THEN
          CALL GetElementEdgeMap( Element, EdgeMap )
          CALL GetElementFaceEdgeMap( Element, FaceEdgeMap ) 
          n = Element % TYPE % NumberOfEdges
       ELSE 
          SELECT CASE( Element % TYPE % ElementCode / 100 )
          CASE(1)
             CYCLE
          CASE(2)
             n = 1
             EdgeMap => LineEdgeMap
             FaceEdgeMap => Null()
          CASE(3)
             n = 3
             EdgeMap => TriEdgeMap
             FaceEdgeMap => Null()
          CASE(4)
             n = 4
             EdgeMap => QuadEdgeMap
             FaceEdgeMap => Null()
          CASE(5)
             n = 6
             EdgeMap => TetraEdgeMap
             FaceEdgeMap => TetraFaceEdgeMap
          CASE(6)
             n = 8
             EdgeMap => PyramidEdgeMap
             FaceEdgeMap => PyramidFaceEdgeMap
          CASE(7)
             n = 9
             EdgeMap => WedgeEdgeMap
             FaceEdgeMap => WedgeFaceEdgeMap
          CASE(8)
             n = 12
             EdgeMap => BrickEdgeMap
             FaceEdgeMap => BrickFaceEdgeMap
          CASE DEFAULT
             CYCLE
             WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
             CALL Fatal('FindMeshEdges',Message)
          END SELECT
       END IF

!      Loop over every edge of every element:
!      --------------------------------------
       DO k=1,n

!         Use MIN(Node1,Node2) as key to hash table:
!         ------------------------------------------
          n1 = Element % NodeIndexes(EdgeMap(k,1))
          n2 = Element % NodeIndexes(EdgeMap(k,2))
          IF ( n1 < n2 ) THEN
             Node1 = n1
             Node2 = n2
          ELSE
             Node1 = n2
             Node2 = n1
          END IF
!
!         Look the edge from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node1 == Node2 ) THEN
                Found = .TRUE.
                Edge = HashPtr % Edge
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO
!
!         Existing edge, update structures:
!         ---------------------------------
          IF ( Found ) THEN
             Element % EdgeIndexes(k) = Edge

             ! Mark edge as an edge of pydamid square face 
             IF (isPPyramid(Element) .AND. k < 5) THEN
                Edges(Edge) % PDefs % pyramidQuadEdge = .TRUE.
             END IF

             IF ( ASSOCIATED(Mesh % Faces) .AND. ASSOCIATED(FaceEdgeMap) ) THEN
               DO ii=1,Element % TYPE % NumberOfFaces
                 Face => Mesh % Faces(Element % FaceIndexes(ii))
                 IF ( .NOT. ASSOCIATED(Face % EdgeIndexes) ) THEN
                   ALLOCATE(Face % EdgeIndexes(Face % TYPE % NumberOfEdges))
                   Face % EdgeIndexes = 0
                 END IF
                 DO jj=1,Face % TYPE % NumberOfEdges
                    IF (FaceEdgeMap(ii,jj) == k) THEN
                       Face % EdgeIndexes(jj) = Edge
                       IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
                         Edges(Edge) % BoundaryInfo % Left => Face
                       ELSE
                         Edges(Edge) % BoundaryInfo % Right => Face
                       END IF
                       EXIT
                    END IF
                 END DO
               END DO
             ELSE
               IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
                 Edges(Edge) % BoundaryInfo % Left  => Element
               ELSE
                 Edges(Edge) % BoundaryInfo % Right => Element
               END IF
             END IF
          ELSE

!            Edge not yet there, create:
!            ---------------------------
             NofEdges = NofEdges + 1
             Edge = NofEdges
             IF( Edge > SIZE( Edges ) ) THEN
               CALL Fatal('FindMeshEdges3D','Number of edges larger than expected!')
             END IF

             Edges(Edge) % ElementIndex = Edge
             Degree = Element % TYPE % BasisFunctionDegree

!            Edge is always a line segment with deg+1 nodes:
!            -----------------------------------------------
             Edges(Edge) % TYPE => GetElementType( 201 + degree, .FALSE.)

             Edges(Edge) % NDOFs  = 0
             IF (Element % NDOFs /= 0) Edges(Edge) % NDOFs = &
                 Element % NDOFs / Element % TYPE % NumberOfNodes * &
                     Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             Edges(Edge) % EdgeIndexes => NULL()
             Edges(Edge) % FaceIndexes => NULL()

             CALL AllocateVector( Edges(Edge) % NodeIndexes, degree + 1 )
             DO n2=1,degree+1
               Edges(Edge) % NodeIndexes(n2) = &
                    Element % NodeIndexes(EdgeMap(k,n2))
             END DO

             Element % EdgeIndexes(k) = Edge
             ALLOCATE( Edges(Edge) % BoundaryInfo )
             Edges(Edge) % BoundaryInfo % Left  => NULL()
             Edges(Edge) % BoundaryInfo % Right => NULL()

             ! Allocate P element definitions 
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
                CALL AllocatePDefinitions(Edges(Edge))
             
                Edges(Edge) % PDefs % P = 0
                Edges(Edge) % PDefs % pyramidQuadEdge = .FALSE.
                ! Here mark edge as edge of pyramid if needed (or set as not)
                IF (isPPyramid(Element) .AND. k < 5) THEN
                   Edges(Edge) % PDefs % pyramidQuadEdge = .TRUE.
                END IF
             ELSE
                NULLIFY( Edges(Edge) % PDefs )
             END IF

             IF ( ASSOCIATED(Mesh % Faces) .AND. ASSOCIATED(FaceEdgeMap) ) THEN
               DO ii=1,Element % TYPE % NumberOfFaces
                 Face => Mesh % Faces(Element % FaceIndexes(ii))
                 IF (.NOT.ASSOCIATED(Face % EdgeIndexes)) THEN
                    ALLOCATE(Face % EdgeIndexes(Face % TYPE % NumberOfEdges))
                    Face % EdgeIndexes = 0
                 END IF
                 DO jj=1,Face % TYPE % NumberOfEdges
                    IF (FaceEdgeMap(ii,jj) == k) THEN
                       Face % EdgeIndexes(jj) = Edge
                       IF (.NOT.ASSOCIATED( Edges(Edge) % BoundaryInfo % Left)) THEN
                         Edges(Edge) % BoundaryInfo % Left => Face
                       ELSE
                         Edges(Edge) % BoundaryInfo % Right => Face
                       END IF
                    END IF
                 END DO
               END DO
             ELSE
                 IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
                   Edges(Edge) % BoundaryInfo % Left  => Element
                 ELSE
                   Edges(Edge) % BoundaryInfo % Right => Element
                 END IF
               END IF

!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr )
             HashPtr % Edge = Edge
             HashPtr % Node1 = Node2
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO
#endif    
    Mesh % NumberOfEdges = NofEdges
    CALL Info('FindMeshEdges3D','Number of edges found: '//TRIM(I2S(NofEdges)),Level=10)
    
!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )
    
    IF (ASSOCIATED(Mesh % Faces)) CALL FixFaceEdges()

    CALL Info('FindMeshEdges3D','All done',Level=20)

CONTAINS 

    SUBROUTINE FixFaceEdges()

      INTEGER :: i,j,k,n,swap,edgeind(4),i1(2),i2(2)

      DO i=1,Mesh % NumberOfFaces
        Face => Mesh % Faces(i)
        n = Face % TYPE % NumberOfEdges
        Edgeind(1:n) = Face % EdgeIndexes(1:n)
        DO j=1,n
          i1 = Mesh % Edges(Edgeind(j)) % NodeIndexes(1:2)
          IF ( i1(1)>i1(2) ) THEN
            swap=i1(1)
            i1(1)=i1(2)
            i1(2)=swap
          END IF
          DO k=1,n
            i2(1) = k
            i2(2) = k+1
            IF ( i2(2)>n ) i2(2)=1
            i2 = Face % NodeIndexes(i2)
            IF ( i2(1)>i2(2) ) THEN
              swap=i2(1)
              i2(1)=i2(2)
              i2(2)=swap
            END IF
            IF ( ALL(i1 == i2) ) THEN
              Face % EdgeIndexes(k) = edgeind(j)
              EXIT
            END IF
          END DO
        END DO
      END DO
    END SUBROUTINE FixFaceEdges
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Finds neighbours of the nodes in given direction.
!> The algorithm finds the neighbour that within 45 degrees of the 
!> given direction has the smallest distance.
!------------------------------------------------------------------------------
  SUBROUTINE FindNeighbourNodes( Mesh,Direction,Neighbours,EndNeighbours)
!------------------------------------------------------------------------------

  TYPE(Mesh_t) , POINTER :: Mesh 
  REAL(KIND=dp) :: Direction(:)
  INTEGER :: Neighbours(:)
  INTEGER, OPTIONAL :: EndNeighbours(:)

  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  REAL(KIND=dp), POINTER :: Distances(:)
  REAL(KIND=dp) :: rn(3), rs(3), ss, sn
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: i,j,k,n,t,DIM,istat

  IF(SIZE(Neighbours) < Mesh % NumberOfNodes) THEN
    CALL Warn('FindNeigbourNodes','SIZE of Neighbours should equal Number of Nodes!')
    RETURN
  END IF


  IF(PRESENT(EndNeighbours)) THEN
    IF(SIZE(EndNeighbours) < Mesh % NumberOfNodes) THEN
      CALL Warn('FindNeigbourNodes','SIZE of EndNeigbours should equal Number of Nodes!')
      RETURN
    END IF
  END IF


  DIM = CoordinateSystemDimension()
  N = Mesh % MaxElementNodes

  CALL AllocateVector( ElementNodes % x, n )
  CALL AllocateVector( ElementNodes % y, n )
  CALL AllocateVector( ElementNodes % z, n )
  CALL AllocateVector( Distances, Mesh % NumberOfNodes )

  Neighbours = 0
  Distances = HUGE(Distances)
 
  rn(1:DIM) = Direction(1:DIM)
  ss = SQRT(SUM(rn(1:DIM)**2))
  rn = rn / ss

  DO t=1,Mesh % NumberOfBulkElements

    CurrentElement => Mesh % Elements(t)
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
  
    ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
    ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
    IF(DIM == 3) THEN
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))
    END IF


    DO i=1,n
      DO j=i+1,n
        rs(1) = ElementNodes % x(j) - ElementNodes % x(i)
        rs(2) = ElementNodes % y(j) - ElementNodes % y(i)
        IF (DIM == 3) THEN
          rs(3) = ElementNodes % z(j) - ElementNodes % z(i)
        END IF
        
        ss = SQRT(SUM(rs(1:DIM)**2))
        sn = SUM(rs(1:DIM)*rn(1:DIM))

        IF(ss < SQRT(2.0) * ABS(sn)) THEN
          IF(sn > 0) THEN
            IF(ss < Distances(NodeIndexes(i))) THEN
              Distances(NodeIndexes(i)) = ss
              Neighbours(NodeIndexes(i)) = NodeIndexes(j)
            END IF
          ELSE
            IF(ss < Distances(NodeIndexes(j))) THEN
              Distances(NodeIndexes(j)) = ss
              Neighbours(NodeIndexes(j)) = NodeIndexes(i)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO

  ! This loop finds the final neighbour in the end of the chain 
  IF(PRESENT(EndNeighbours)) THEN
    EndNeighbours = Neighbours

    DO t=1,Mesh%NumberOfNodes
      j = Neighbours(t)
      DO WHILE(j /= 0)
        EndNeighbours(t) = j
        j = Neighbours(j)
      END DO
    END DO
  END IF
  DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z, Distances)
!------------------------------------------------------------------------------
END SUBROUTINE FindNeighbourNodes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE UpdateSolverMesh( Solver, Mesh )
!------------------------------------------------------------------------------
     TYPE( Mesh_t ), POINTER :: Mesh
     TYPE( Solver_t ), TARGET :: Solver
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,n1,n2,DOFs
     LOGICAL :: Found, OptimizeBandwidth
     TYPE(Matrix_t), POINTER   :: Matrix
     REAL(KIND=dp), POINTER :: Work(:)
     INTEGER, POINTER :: Permutation(:)
     TYPE(Variable_t), POINTER :: TimeVar, SaveVar, Var
     CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
     SaveVar => Solver % Variable
     DOFs = SaveVar % DOFs

     Solver % Mesh => Mesh
     CALL SetCurrentMesh( CurrentModel, Mesh )
!
!    Create matrix and variable structures for
!    current equation on the new mesh:
!    -----------------------------------------
     Solver % Variable => VariableGet( Mesh % Variables, &
        Solver % Variable % Name, ThisOnly = .FALSE. )

     CALL AllocateVector( Permutation, SIZE(Solver % Variable % Perm) )

     OptimizeBandwidth = ListGetLogical( Solver % Values, 'Optimize Bandwidth', Found )
     IF ( .NOT. Found ) OptimizeBandwidth = .TRUE.

     Matrix => CreateMatrix( CurrentModel, Solver, &
        Mesh, Permutation, DOFs, MATRIX_CRS, OptimizeBandwidth, &
        ListGetString( Solver % Values, 'Equation' ) )

     Matrix % Symmetric = ListGetLogical( Solver % Values, &
             'Linear System Symmetric', Found )

     Matrix % Lumped = ListGetLogical( Solver % Values, &
             'Lumped Mass Matrix', Found )

     ALLOCATE( Work(SIZE(Solver % Variable % Values)) )
     Work = Solver % Variable % Values
     DO k=0,DOFs-1
        DO i=1,SIZE(Permutation)
           IF ( Permutation(i) > 0 ) THEN
              Solver % Variable % Values( DOFs*Permutation(i)-k ) = &
                 Work( DOFs*Solver % Variable % Perm(i)-k )
           END IF
        END DO
     END DO

     IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        DO j=1,SIZE(Solver % Variable % PrevValues,2)
           Work = Solver % Variable % PrevValues(:,j)
           DO k=0,DOFs-1
              DO i=1,SIZE(Permutation)
                 IF ( Permutation(i) > 0 ) THEN
                    Solver % Variable % PrevValues( DOFs*Permutation(i) - k,j ) =  &
                        Work( DOFs * Solver % Variable % Perm(i) - k )
                  END IF
              END DO
           END DO
        END DO
     END IF
     DEALLOCATE( Work )

     Solver % Variable % Perm = Permutation
     Solver % Variable % Solver => Solver

     DEALLOCATE( Permutation )
     CALL AllocateVector( Matrix % RHS, Matrix % NumberOfRows )

     IF ( ASSOCIATED(SaveVar % EigenValues) ) THEN
        n = SIZE(SaveVar % EigenValues)

        IF ( n > 0 ) THEN
           Solver % NOFEigenValues = n
           CALL AllocateVector( Solver % Variable % EigenValues,n )
           CALL AllocateArray( Solver % Variable % EigenVectors, n, &
                    SIZE(Solver % Variable % Values) ) 

           IF( Solver % Variable % Dofs > 1 ) THEN
             DO k=1,Solver % Variable % DOFs
               str = ComponentName( Solver % Variable % Name, k )
               Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )
               IF ( ASSOCIATED( Var ) ) THEN
                 Var % EigenValues => Solver % Variable % EigenValues
                 Var % EigenVectors =>  & 
                     Solver % Variable % EigenVectors(:,k::Solver % Variable % DOFs )
               END IF
             END DO
           END IF
           
           Solver % Variable % EigenValues  = 0.0d0
           Solver % Variable % EigenVectors = 0.0d0

           CALL AllocateVector( Matrix % MassValues, SIZE(Matrix % Values) )
           Matrix % MassValues = 0.0d0
        END IF
     ELSE IF ( ASSOCIATED( Solver % Matrix ) ) THEN
        IF( ASSOCIATED( Solver % Matrix % Force) ) THEN
           n1 = Matrix % NumberOFRows
           n2 = SIZE(Solver % Matrix % Force,2)
           ALLOCATE(Matrix % Force(n1,n2))
           Matrix % Force = 0.0d0
        END IF
     END IF

     Solver % Matrix => Matrix
     Solver % Mesh % Changed = .TRUE.

!------------------------------------------------------------------------------
  END SUBROUTINE UpdateSolverMesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Split a mesh equally to smaller pieces by performing a uniform split.
!> Also known as mesh multiplication. A 2D element splits into 4 elements of
!> same form, and 3D element into 8 elements. 
!> Currently works only for linear elements.
!------------------------------------------------------------------------------
  FUNCTION SplitMeshEqual(Mesh,h) RESULT( NewMesh )
!------------------------------------------------------------------------------
    REAL(KIND=dp), OPTIONAL :: h(:)
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: u(:),v(:),w(:),x(:),y(:),z(:),xh(:)
    INTEGER :: i, j, k, n, NewElCnt, NodeCnt, EdgeCnt, FaceCnt, Node, ParentId, Diag, NodeIt
    LOGICAL :: Found, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Eptr,Eparent,Face,Faces(:)
    INTEGER, POINTER :: Child(:,:)
    INTEGER :: n1,n2,n3,EoldNodes(4),FaceNodes(4),EdgeNodes(2) ! Only linears so far
    INTEGER :: FaceNumber,Edge1,Edge2,Edge3,Edge4,Node12,Node23,Node34,Node41,Node31
    REAL(KIND=dp) :: dxyz(3,3),Dist(3),r,s,t,h1,h2
    TYPE(PElementDefs_t), POINTER :: PDefs
    INTEGER :: ierr, ParTmp(6), ParSizes(6)
    INTEGER, ALLOCATABLE :: FacePerm(:), BulkPerm(:)
    LOGICAL :: Parallel
!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Mesh ) ) RETURN

    CALL Info( 'SplitMeshEqual', 'Mesh splitting works for first order elements 303, 404, 504, (706) and 808.', Level = 6 )

    DO i=1,Mesh % NumberOfBulkElements
      SELECT CASE(Mesh % Elements(i) % TYPE % ElementCode/100)
      CASE(6)
        CALL Fatal('SplitMeshEqual','Pyramids not supported, sorry.')
      END SELECT
    END DO

    NewMesh => AllocateMesh()

    NewMesh % SingleMesh = Mesh % SingleMesh
    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. NewMesh % SingleMesh )

    
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT.EdgesPresent) CALL FindMeshEdges( Mesh )

    CALL ResetTimer('SplitMeshEqual')

    CALL Info( 'SplitMeshEqual', '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( 'SplitMeshEqual', Message, Level=6 )
!
!   Update nodal coordinates:
!   -------------------------
    NodeCnt = Mesh % NumberOfNodes + Mesh % NumberOfEdges
!
!   For quad faces add one node in the center:
!   ------------------------
    ALLOCATE(FacePerm(Mesh % NumberOfFaces)); FacePerm = 0
    FaceCnt = 0
    DO i = 1, Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       IF( Face % TYPE % NumberOfNodes == 4 ) THEN
         NodeCnt = NodeCnt+1
         FaceCnt = FaceCnt+1
         FacePerm(i) = NodeCnt
       END IF
    END DO
    
    WRITE( Message, * ) 'Added nodes in the center of faces : ', FaceCnt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )
!
!   For quads and bricks, count centerpoints:
!   -----------------------------------------
    NodeIt = 0
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(4,8)
          NodeCnt = NodeCnt + 1
          NodeIt = NodeIt + 1
       END SELECT
    END DO
    
    WRITE( Message, * ) 'Added nodes in the center of bulks : ', NodeIt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )
!
!   new mesh nodecoordinate arrays:
!   -------------------------------
    CALL AllocateVector( NewMesh % Nodes % x, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % y, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % z, NodeCnt )

!   shortcuts (u,v,w) old mesh  nodes,
!   (x,y,z) new mesh nodes:
!   ----------------------------------
    u => Mesh % Nodes % x
    v => Mesh % Nodes % y
    w => Mesh % Nodes % z

    x => NewMesh % Nodes % x
    y => NewMesh % Nodes % y
    z => NewMesh % Nodes % z
!
!   new mesh includes old mesh nodes:
!   ----------------------------------
    x(1:Mesh % NumberOfNodes) = u
    y(1:Mesh % NumberOfNodes) = v
    z(1:Mesh % NumberOfNodes) = w

! what is h? - pointer to nodal element size
    IF (PRESENT(h)) THEN
      ALLOCATE(xh(SIZE(x)))
      xh(1:SIZE(h)) = h
    END IF
!
!   add edge centers:
!   -----------------
    j =  Mesh % NumberOfNodes
    DO i=1,Mesh % NumberOfEdges
       j = j + 1
       Edge => Mesh % Edges(i)
       k = Edge % TYPE % NumberOfNodes
       IF (PRESENT(h)) THEN
         h1=h(Edge % NodeIndexes(1))
         h2=h(Edge % NodeIndexes(2))
         r=1._dp/(1+h1/h2)
         x(j) = r*u(Edge%NodeIndexes(1))+(1-r)*u(Edge%NodeIndexes(2))
         y(j) = r*v(Edge%NodeIndexes(1))+(1-r)*v(Edge%NodeIndexes(2))
         z(j) = r*w(Edge%NodeIndexes(1))+(1-r)*w(Edge%NodeIndexes(2))
         xh(j)=r*h1+(1-r)*h2
       ELSE
         x(j) = SUM(u(Edge % NodeIndexes))/k
         y(j) = SUM(v(Edge % NodeIndexes))/k
         z(j) = SUM(w(Edge % NodeIndexes))/k
       END IF
    END DO
    
    CALL Info('SplitMeshEqual','Added edge centers to the nodes list.', Level=10 )  
!
!   add quad face centers for bricks and prisms(wedges):
!   ----------------------------
    j = Mesh % NumberOfNodes + Mesh % NumberOfEdges
    DO i=1,Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       k = Face % TYPE % NumberOfNodes
       IF( k == 4 ) THEN
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Face % EdgeIndexes(2))
            h2=xh(n+Face % EdgeIndexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Face % EdgeIndexes(3))
            h2=xh(n+Face % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Face,u(Face % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Face,v(Face % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Face,w(Face % NodeIndexes),r,s)
            xh(j) = InterpolateInElement2D(Face,h(Face % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Face % NodeIndexes))/k
            y(j) = SUM(v(Face % NodeIndexes))/k
            z(j) = SUM(w(Face % NodeIndexes))/k
          END IF
       END IF
    END DO
    
    CALL Info('SplitMeshEqual','Added face centers to the nodes list.', Level=10 )
!
!   add centerpoint for quads & bricks:
!   -----------------------------------
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       k = Eold % TYPE % NumberOfNodes
       SELECT CASE( Eold % TYPE % ElementCode / 100 )

       CASE(4)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Eold % Edgeindexes(2))
            h2=xh(n+Eold % Edgeindexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Eold % EdgeIndexes(3))
            h2=xh(n+Eold % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Eold,u(Eold % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Eold,v(Eold % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Eold,w(Eold % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       CASE(8)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes+Mesh % NumberOfEdges
            h1=xh(n+Eold % FaceIndexes(4))
            h2=xh(n+Eold % FaceIndexes(6))
            r=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(5))
            h2=xh(n+Eold % FaceIndexes(3))
            s=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(2))
            h2=xh(n+Eold % FaceIndexes(1))
            t=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement3D(Eold,u(Eold % NodeIndexes),r,s,t)
            y(j) = InterpolateInElement3D(Eold,v(Eold % NodeIndexes),r,s,t)
            z(j) = InterpolateInElement3D(Eold,w(Eold % NodeIndexes),r,s,t)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       END SELECT
    END DO
!
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MinEdgeDOFs = Mesh % MinEdgeDOFs
    NewMesh % MinFaceDOFs = Mesh % MinFaceDOFs
    NewMesh % MaxEdgeDOFs = Mesh % MaxEdgeDOFs
    NewMesh % MaxFaceDOFs = Mesh % MaxFaceDOFs
    NewMesh % MaxElementDOFs = Mesh % MaxElementDOFs
    NewMesh % MeshDim = Mesh % MeshDim

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt
!
!   Update bulk elements:
!   =====================
!
!   First count new elements:
!   -------------------------
    NewElCnt = 0
    DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode/100 )

!      Each element will be divided into 2**Dim new elements:
!      ------------------------------------------------------
       CASE(2)
          NewElCnt = NewElCnt + 2 ! lines
       CASE(3)
          NewElCnt = NewElCnt + 4 ! trias
       CASE(4)
          NewElCnt = NewElCnt + 4 ! quads
       CASE(5)
          NewElCnt = NewElCnt + 8 ! tetras
       CASE(7)
          NewElCnt = NewElCnt + 8 ! prisms (wedges)
       CASE(8)
          NewElCnt = NewElCnt + 8 ! hexas
       END SELECT
    END DO

    WRITE( Message, * ) 'Count of new elements : ', NewElCnt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info('SplitMeshEqual','New mesh allocated.', Level=10 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 8 )
    CALL Info('SplitMeshEqual','Array for bulk elements allocated.', Level=10 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes
    EdgeCnt = Mesh % NumberOfEdges

!
!   Index to old quad/hexa centerpoint node in the new mesh nodal arrays:
!   ---------------------------------------------------------------------
    Node = NodeCnt + EdgeCnt + FaceCnt
!
!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)

       SELECT CASE( Eold % TYPE % ElementCode )
       CASE(303)
!
!         Split triangle to four triangles from
!         edge centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,2) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,3) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,4) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt

       CASE(404)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split quad to four new quads from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)


       CASE(504)
!
!         Split tetra to 8 new elements from
!         corners and edge centerpoints:
!         ----------------------------------
!
!         1st new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt

!         Then the annoying part; we still have to split the
!         remaining octahedron into four elements. This can
!         be done in three ways of which only one preserves
!         the minimum angle condition (Delaunay splitting):
!         --------------------------------------------------
          dxyz(1,1) = x(Eold % EdgeIndexes(4) + NodeCnt) &
                    - x(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(2,1) = y(Eold % EdgeIndexes(4) + NodeCnt) &
                    - y(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(3,1) = z(Eold % EdgeIndexes(4) + NodeCnt) &
                    - z(Eold % EdgeIndexes(2) + NodeCnt)

          dxyz(1,2) = x(Eold % EdgeIndexes(5) + NodeCnt) &
                    - x(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(2,2) = y(Eold % EdgeIndexes(5) + NodeCnt) &
                    - y(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(3,2) = z(Eold % EdgeIndexes(5) + NodeCnt) &
                    - z(Eold % EdgeIndexes(3) + NodeCnt)

          dxyz(1,3) = x(Eold % EdgeIndexes(6) + NodeCnt) &
                    - x(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(2,3) = y(Eold % EdgeIndexes(6) + NodeCnt) &
                    - y(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(3,3) = z(Eold % EdgeIndexes(6) + NodeCnt) &
                    - z(Eold % EdgeIndexes(1) + NodeCnt)

          Dist(1) = SQRT( dxyz(1,1)**2 + dxyz(2,1)**2 + dxyz(3,1)**2 )
          Dist(2) = SQRT( dxyz(1,2)**2 + dxyz(2,2)**2 + dxyz(3,2)**2 )
          Dist(3) = SQRT( dxyz(1,3)**2 + dxyz(2,3)**2 + dxyz(3,3)**2 )

          Diag = 1  ! The default diagonal for splitting is between edges 2-4
          IF (Dist(2) < Dist(1) .AND. Dist(2) < Dist(3)) Diag = 2 ! Edges 3-5
          IF (Dist(3) < Dist(1) .AND. Dist(3) < Dist(2)) Diag = 3 ! Edges 1-6

          SELECT CASE( Diag )
          CASE(1)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
          CASE(2)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
          CASE(3)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt

          END SELECT


       CASE(706)
!
!         Split prism to 8 new prism from edge
!         centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt 
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt 
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt 
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))

!
!         3rd new element (near node 3)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(9) + NodeCnt

!
!         4th new element (bottom center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt

!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt

!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
!
!         8th new element (top half, center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt



       CASE(808)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split brick to 8 new bricks from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 8)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(7) = Node
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(6))
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(8) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(6) = Node
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(12)+ NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(5) = Node
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(5))
!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(8) + NodeCnt
!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Node
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(2))
!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(12)+ NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(8) = Eold % NodeIndexes(8)
!
!         8th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(7) = Eold % NodeIndexes(7)
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(7) + NodeCnt

       CASE DEFAULT
          WRITE( Message,* ) 'Element type ', Eold % TYPE % ElementCode, &
              ' not supported by the multigrid solver.'
          CALL Fatal( 'SplitMeshEqual', Message )
       END SELECT
    END DO

!
!   Update new mesh element counts:
!   -------------------------------
    NewMesh % NumberOfBulkElements = NewElCnt

!
!   Update boundary elements:
!   NOTE: Internal boundaries not taken care of...:!!!!
!   ---------------------------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements

       j = i + Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(j)
!
!      get parent of the boundary element:
!      -----------------------------------
       Eparent => Eold % BoundaryInfo % Left
       IF ( .NOT.ASSOCIATED(Eparent) ) &
          eParent => Eold % BoundaryInfo % Right
       IF ( .NOT. ASSOCIATED( Eparent ) ) CYCLE

       ParentId = Eparent % ElementIndex

       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(2)
!
!         Line segments:
!         ==============
!
!         which edge of the parent element are we ?
!         -----------------------------------------
          DO Edge1=1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             IF ( Eold % NodeIndexes(1) == Edge % NodeIndexes(1) .AND. &
                  Eold % NodeIndexes(2) == Edge % NodeIndexes(2) .OR.  &
                  Eold % NodeIndexes(2) == Edge % NodeIndexes(1) .AND. &
                  Eold % NodeIndexes(1) == Edge % NodeIndexes(2) ) EXIT
          END DO
!
!         index of the old edge centerpoint in the
!         new mesh nodal arrays:
!         ----------------------------------------
          Node = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,4
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found = .FALSE.
             DO k=1,n-1
                IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k+1) .OR.  &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(1) == Eptr % NodeIndexes(k+1) ) THEN
                   Found = .TRUE.
                   EXIT
                END IF
             END DO
             IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(1) .OR.  &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(1) == Eptr % NodeIndexes(1) ) THEN
                Found = .TRUE.
             END IF
             IF ( Found ) EXIT
          END DO
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,4
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found = .FALSE.
             DO k=1,n-1
                IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k+1) .OR.  &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(1) == Eptr % NodeIndexes(k+1) ) THEN
                   Found = .TRUE.
                   EXIT
                END IF
             END DO
             IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(1) .OR.  &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(1) == Eptr % NodeIndexes(1) ) THEN
                Found = .TRUE.
             END IF
             IF ( Found ) EXIT
          END DO
          Enew % BoundaryInfo % Left => Eptr

       CASE(3)
!
!         Trias:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:3) = Eold % NodeIndexes(1:3)
          CALL sort( 3, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:3) = Face % NodeIndexes(1:3)
             CALL sort( 3, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) ) EXIT

          END DO
!
!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node31 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node31
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr

       CASE(4)
!
!         Quads:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:4) = Eold % NodeIndexes(1:4)
          CALL sort( 4, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:4) = Face % NodeIndexes(1:4)
             CALL sort( 4, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) .AND. &
                  EoldNodes(4) == FaceNodes(4) ) EXIT

          END DO

!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Fourth edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          DO Edge4 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge4) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node = FacePerm(Eparent % FaceIndexes(FaceNumber)) ! faces mid-point
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node34 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
          Node41 = Eparent % EdgeIndexes(Edge4) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Node41
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 )  CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          Enew % NodeIndexes(4) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node41
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Node34
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Node34
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
       END SELECT
    END DO

!
!   Update new mesh boundary element counts:
!   ----------------------------------------
    NewMesh % NumberOfBoundaryElements = NewElCnt - &
            NewMesh % NumberOfBulkElements
    NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs
    NewMesh % MaxElementNodes = Mesh % MaxElementNodes

    j = 0
    DO i=1,NewMesh % NumberOfBulkElements+NewMesh % NumberOfBoundaryElements
      Enew => NewMesh % Elements(i)

      IF ( Enew % DGDOFs>0 ) THEN
        ALLOCATE(Enew % DGIndexes(Enew % DGDOFs))
        DO k=1,Enew % DGDOFs
          j = j + 1
          Enew % DGIndexes(k)=j
        END DO
      ELSE
        Enew % DGIndexes=>NULL()
      END IF

      IF (i<=NewMesh % NumberOfBulkElements) THEN
         PDefs => Enew % PDefs

         IF(ASSOCIATED(PDefs)) THEN
           CALL AllocatePDefinitions(Enew)
           Enew % PDefs = PDefs

           ! All elements in actual mesh are not edges
           Enew % PDefs % pyramidQuadEdge = .FALSE.
           Enew % PDefs % isEdge = .FALSE.

           ! If element is of type tetrahedron and is a p element,
           ! do the Ainsworth & Coyle trick
           IF (Enew % TYPE % ElementCode == 504) CALL ConvertToACTetra(Enew)
            CALL GetRefPElementNodes( Enew % Type,  Enew % Type % NodeU, &
                 Enew % Type % NodeV, Enew % Type % NodeW )
         END IF
      ELSE
        Enew % PDefs=>NULL()
      END IF
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()
      Enew % BubbleIndexes => NULL()
    END DO

    CALL Info( 'SplitMeshEqual', '******** New mesh ********', Level=6 )
    WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
    CALL Info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
    CALL Info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
    CALL Info( 'SplitMeshEqual', Message, Level=6 )


    ! Information of the new system size, also in parallel
    !----------------------------------------------------------------------
    ParTmp(1) = Mesh % NumberOfNodes
    ParTmp(2) = Mesh % NumberOfBulkElements
    ParTmp(3) = Mesh % NumberOfBoundaryElements
    ParTmp(4) = NewMesh % NumberOfNodes
    ParTmp(5) = NewMesh % NumberOfBulkElements
    ParTmp(6) = NewMesh % NumberOfBoundaryElements

    IF( .FALSE. .AND. Parallel ) THEN
      CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)

      CALL Info('SplitMeshEqual','Information on parallel mesh sizes')
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(1),' nodes'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(2),' bulk elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(3),' boundary elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(4),' nodes'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(5),' bulk elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(6),' boundary elements'
      CALL Info('SplitMeshEqual',Message)
    END IF


    CALL CheckTimer('SplitMeshEqual',Delete=.TRUE.)

!
!   Update structures needed for parallel execution:
!   ------------------------------------------------
    IF( Parallel ) THEN
      CALL UpdateParallelMesh( Mesh, NewMesh )
    END IF
!
!
!   Finalize:
!   ---------
    DEALLOCATE( Child )
    IF(.NOT.EdgesPresent) THEN
      CALL Info('SplitMeshEqual','Releasing edges from the old mesh as they are not needed!',Level=20)
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    ELSE
      CALL Info('SplitMeshEqual','Generating edges in the new mesh as thet were present in the old!',Level=20)
      CALL FindMeshEdges( NewMesh )
    END IF

    ! Our boundary may be a circle, cylider or sphere surface.
    ! Honor those shapes when splitting the mesh!
    CALL FollowCurvedBoundary( CurrentModel, NewMesh, .FALSE. ) 

    
!call writemeshtodisk( NewMesh, "." )
!stop
CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE UpdateParallelMesh( Mesh, NewMesh )
!------------------------------------------------------------------------------
       TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
       TYPE(Element_t), POINTER :: Edge, Face, Element, BoundaryElement
       INTEGER :: i,j,k,l,m,n,p,q, istat
       INTEGER, POINTER :: IntCnts(:),IntArray(:),Reorder(:)
       INTEGER, ALLOCATABLE :: list1(:), list2(:)
       LOGICAL, ALLOCATABLE :: InterfaceTag(:)

       INTEGER :: jedges
       LOGICAL :: Found
!------------------------------------------------------------------------------

!
!      Update mesh interfaces for parallel execution.
!      ==============================================
!
!      Try to get an agreement about the  global numbering
!      of new mesh nodes among set of processes solving
!      this specific eq. Also allocate and generate
!      all other control information needed in parallel
!      execution:
!      ----------------------------------------------------
       n = NewMesh % NumberOfNodes
       ALLOCATE( NewMesh % ParallelInfo % NeighbourList(n), stat=istat )
       IF ( istat /= 0 ) &
         CALL Fatal( 'UpdateParallelMesh', 'Allocate error.' )
       CALL AllocateVector( NewMesh % ParallelInfo % NodeInterface,n  )
       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )

       DO i=1,n
          NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % NodeInterface = .FALSE.
       NewMesh % ParallelInfo % NodeInterface(1:n) = Mesh % ParallelInfo % NodeInterface

       NewMesh % ParallelInfo % GlobalDOFs = 0
       NewMesh % ParallelInfo % GlobalDOFs(1:n) = &
          Mesh % ParallelInfo % GlobalDOFs
!
!      My theory is, that a new node will be an
!      interface node only if all the edge or face
!      nodes which contribute to its existence are
!      interface nodes (the code immediately below
!      will only count sizes):
!      -------------------------------------------
!

       ! New version based on edges and faces (2. March 2007):
       !=====================================================
       SELECT CASE( CoordinateSystemDimension() )
          
       CASE(2)
          !
          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % NodeInterface(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'
          !
          ! Determine possible interface edges:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfEdges ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfEdges
             Edge => Mesh % Edges(i)
             IF( ASSOCIATED(Edge % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Edge % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % NodeInterface( Edge % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          !
          ! Eliminate false positives based on BoundaryElement -data:
          !----------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements( Mesh % NumberOfBulkElements + i )
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED( Element ) ) &
                  Element => BoundaryElement % BoundaryInfo % Right
             IF( .NOT.ASSOCIATED( Element ) ) CYCLE
             IF( .NOT.ASSOCIATED( Element % EdgeIndexes ) ) CYCLE
             
             ALLOCATE( list1( SIZE( BoundaryElement % NodeIndexes )))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort( SIZE(list1), list1 )
             
             DO j = 1,Element % TYPE % NumberOfEdges
                k = Element % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                IF( SIZE( Edge % NodeIndexes ) /= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Edge % NodeIndexes )))
                list2 = Edge % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO

                DEALLOCATE(list2)
                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfEdges
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Edge => Mesh % Edges(i)
             
             ! This is just for the edge count:
             !---------------------------------
             IF( NewMesh % ParallelInfo % NodeInterface( Mesh % NumberOfNodes + i) ) CYCLE
             
             ! Mark interface nodes and count edges:
             !--------------------------------------
             NewMesh % ParallelInfo % NodeInterface( Mesh % NumberOfNodes + i) = .TRUE.
             p = p+1

          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 2*p ! check
          
       CASE(3)

          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % NodeInterface(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'

          ! Determine possible interface faces:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfFaces ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( ASSOCIATED(Face % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Face % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % NodeInterface( Face % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          
          ! Eliminate false interface faces based on BoundaryElement -data:
          !----------------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements(Mesh % NumberOfBulkElements+i)
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED(Element) ) &
                Element => BoundaryElement % BoundaryInfo % Right
              IF( .NOT.ASSOCIATED(Element) ) CYCLE
              IF( .NOT.ASSOCIATED(Element % FaceIndexes) ) CYCLE
             
             ALLOCATE(list1(SIZE(BoundaryElement % NodeIndexes)))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort(SIZE(list1),list1)
             
             DO j = 1,Element % TYPE % NumberOfFaces
                k = Element % FaceIndexes(j)
                Face => Mesh % Faces(k)
                IF(SIZE(Face % NodeIndexes)/= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Face % NodeIndexes )))
                list2 = Face % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO
                
                DEALLOCATE(list2)

                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Count interface faces:
          !-----------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( InterfaceTag(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface faces'
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Face => Mesh % Faces(i)
             
             DO j = 1,SIZE( Face % EdgeIndexes )
                k = Face % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                
                ! This is just for the edge count:
                !---------------------------------
                IF( NewMesh % ParallelInfo % NodeInterface( Mesh % NumberOfNodes + k) ) CYCLE
                
                ! Mark interface nodes and count edges:
                !--------------------------------------
                NewMesh % ParallelInfo % NodeInterface( Mesh % NumberOfNodes + k) = .TRUE.
                p = p+1
             END DO
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 3*p ! check
          
       END SELECT

!======================================================================================================
       j = p
       jedges = p

!      For bricks, check also the faces:
!      ---------------------------------
       DO i = 1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i) 
          IF( Face % TYPE % NumberOfNodes == 4 ) THEN
             IF ( ALL( Mesh % ParallelInfo % NodeInterface( Face % NodeIndexes ) ) ) THEN
                NewMesh % ParallelInfo % NodeInterface( Mesh % NumberOfNodes &
                     + Mesh % NumberOfEdges + i ) = .TRUE.
                j = j + 1
                k = k + Face % TYPE % NumberOfNodes
             END IF
          END IF
       END DO

!      CALL AllocateVector( IntCnts,  j )
!      CALL AllocateVector( IntArray, k )
!
!      Old mesh nodes were copied as is...
!
       IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % Neighbourlist ) ) THEN
         CALL Fatal('UpdateParallelMesh','Original mesh has no NeighbourList!')
       END IF
       
       DO i=1,Mesh % NumberOfNodes
          CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours, &
                SIZE( Mesh % ParallelInfo % Neighbourlist(i) % Neighbours) )

          NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
             Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       END DO
!
!      Take care of the new mesh internal nodes.
!      Parallel global numbering will take care
!      of the interface nodes:
!      ----------------------------------------
       DO i=Mesh % NumberOfNodes+1, NewMesh % NumberOfNodes
          IF ( .NOT. NewMesh % ParallelInfo % NodeInterface(i) ) THEN
            CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours,1 )
            NewMesh % ParallelInfo % NeighbourList(i) %  Neighbours(1) = ParEnv % MyPE
          END IF
       END DO
!
!      Copy global indices of edge and/or face nodes
!      to temporary work arrays:
!      ---------------------------------------------
!
! check also this:
!      j = 0
!      k = 0
!      DO i = 1,Mesh % NumberOfEdges
!         Edge => Mesh % Edges(i)
!         
!         ! Added check for parent elements 25.2.2007:
!         Found = .NOT.( ASSOCIATED(edge % boundaryinfo % left) &
!              .AND.  ASSOCIATED(edge % boundaryinfo % right) )
!         
!         IF ( ALL(Mesh % ParallelInfo % NodeInterface(Edge % NodeIndexes)) .AND. Found ) THEN
!            j = j + 1
!            IntCnts(j) = Edge % TYPE % NumberOfNodes
!            IntArray( k+1:k+IntCnts(j) ) = &
!                 Mesh % Parallelinfo % GlobalDOFs(Edge % NodeIndexes)
!            CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!            k = k + IntCnts(j)
!         END IF
!      END DO
!      !
!      ! For bricks, check also the faces:
!      ! ---------------------------------
!      DO i = 1,Mesh % NumberOfFaces
!         Face => Mesh % Faces(i)
!         IF( Face % TYPE % NumberOfNodes == 4 ) THEN
!            IF ( ALL( Mesh % ParallelInfo % NodeInterface(Face % NodeIndexes) ) ) THEN
!               j = j + 1
!               IntCnts(j) = Face % TYPE % NumberOfNodes
!               IntArray(k+1:k+IntCnts(j)) = &
!                    Mesh % ParallelInfo % GlobalDOFs(Face % NodeIndexes)
!               CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!               k = k + IntCnts(j)
!            END IF
!         END IF
!      END DO
!
!      Finally the beef, do the exchange of new
!      interfaces. The parallel global numbering
!      subroutine will also do reordering of the
!      nodes, hence the reorder array:
!      -------------------------------------------
       CALL AllocateVector( Reorder, NewMesh % NumberOfNodes )
       Reorder = [ (i, i=1,NewMesh % NumberOfNodes) ]

       k = NewMesh % Nodes % NumberOfNodes - Mesh % Nodes % NumberOfNodes
       CALL ParallelGlobalNumbering( NewMesh, Mesh, k, Reorder )

!      Account for the reordering of the nodes:
!      ----------------------------------------
       DO i=1,NewMesh % NumberOfBulkElements + &
            NewMesh % NumberOfBoundaryElements
          NewMesh % Elements(i) % NodeIndexes = &
              Reorder( NewMesh % Elements(i) % NodeIndexes )
       END DO

!      DEALLOCATE( IntCnts, IntArray, Reorder )
!      DEALLOCATE( Reorder )
!------------------------------------------------------------------------------
    END SUBROUTINE UpdateParallelMesh
  END FUNCTION SplitMeshEqual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SetCurrentMesh( Model, Mesh )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t),  POINTER :: Mesh
!------------------------------------------------------------------------------
    Model % Variables => Mesh % Variables

    Model % Mesh  => Mesh
    Model % Nodes => Mesh % Nodes
    Model % NumberOfNodes = Mesh % NumberOfNodes
    Model % Nodes % NumberOfNodes = Mesh % NumberOfNodes

    Model % Elements => Mesh % Elements
    Model % MaxElementNodes = Mesh % MaxElementNodes
    Model % NumberOfBulkElements = Mesh % NumberOfBulkElements
    Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
  END SUBROUTINE SetCurrentMesh
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------
  SUBROUTINE DisplaceMesh( Mesh, Update, sgn, Perm, DOFs, StabRecomp, UpdateDirs )
!----------------------------------------------------------------------------------
    TYPE(Mesh_t) , POINTER :: Mesh 
    REAL(KIND=dp) :: Update(:)
    INTEGER :: DOFs,sgn,Perm(:)
    LOGICAL, OPTIONAL :: StabRecomp
    INTEGER, OPTIONAL :: UpdateDirs

    INTEGER :: i,k,dim
    LOGICAL :: StabFlag

    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element

    IF ( PRESENT( UpdateDirs ) ) THEN
      dim = UpdateDirs
    ELSE
      dim = DOFs
    END IF

    DO i=1,MIN( SIZE(Perm), SIZE(Mesh % Nodes % x) )
       k = Perm(i)
       IF ( k > 0 ) THEN
         k = DOFs * (k-1)
         Mesh % Nodes % x(i)   = Mesh % Nodes % x(i) + sgn * Update(k+1)
         IF ( dim > 1 ) &
           Mesh % Nodes % y(i) = Mesh % Nodes % y(i) + sgn * Update(k+2)
         IF ( dim > 2 ) &
           Mesh % Nodes % z(i) = Mesh % Nodes % z(i) + sgn * Update(k+3)
        END IF
    END DO

    StabFlag = .TRUE.
    IF ( PRESENT( StabRecomp ) ) StabFlag = StabRecomp

    IF ( sgn == 1 .AND. StabFlag ) THEN
       k = Mesh % MaxElementDOFs
       CALL AllocateVector( ElementNodes % x,k )
       CALL AllocateVector( ElementNodes % y,k )
       CALL AllocateVector( ElementNodes % z,k )

       DO i=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(i)
          IF ( ANY( Perm( Element % NodeIndexes ) == 0 ) ) CYCLE

          k = Element % TYPE % NumberOfNodes
          ElementNodes % x(1:k) = Mesh % Nodes % x(Element % NodeIndexes)
          ElementNodes % y(1:k) = Mesh % Nodes % y(Element % NodeIndexes)
          ElementNodes % z(1:k) = Mesh % Nodes % z(Element % NodeIndexes)
          IF ( Mesh % Stabilize ) THEN
             CALL StabParam( Element,ElementNodes,k, &
                          Element % StabilizationMk, Element % Hk )
          ELSE
             Element % hK = ElementDiameter( Element, ElementNodes )
          END IF
       END DO

       DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE DisplaceMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Convert tetrahedral element to Ainsworth & Coyle type tetrahedron.
!------------------------------------------------------------------------------
  SUBROUTINE ConvertToACTetra( Tetra )
!------------------------------------------------------------------------------
    USE PElementMaps, ONLY : getTetraEdgeMap, getTetraFaceMap
    IMPLICIT NONE
    
    TYPE(Element_t), POINTER :: Tetra  !< Tetrahedral element to convert
!------------------------------------------------------------------------------
    INTEGER :: i, globalMin, globalMax, globalMinI
    INTEGER, DIMENSION(3) :: face, globalFace
    INTRINSIC MIN, MAX, CSHIFT

    ! Sanity check
    IF (Tetra % TYPE % ElementCode /= 504 .OR. &
         .NOT. ASSOCIATED(Tetra % PDefs)) THEN
       CALL Warn('MeshUtils::ConvertToACTetra','Element to convert not p tetrahedron!')
       RETURN
    END IF    
   
    ! Find global min and max vertices
    globalMin = Tetra % NodeIndexes(1)
    globalMinI = 1
    globalMax = Tetra % NodeIndexes(1)
    DO i=2,4
       ! Find min
       IF (globalMin > Tetra % NodeIndexes(i)) THEN
          globalMin = Tetra % NodeIndexes(i)
          globalMinI = i
       ELSE IF (globalMax < Tetra % NodeIndexes(i)) THEN
          globalMax = Tetra % NodeIndexes(i)
       END IF
    END DO
    
    ! Get face containing global min (either face 1 or 2)
    IF (globalMinI == 4) THEN
       face = getTetraFaceMap(2)
    ELSE
       face = getTetraFaceMap(1)
    END IF
    globalFace(1:3) = Tetra % NodeIndexes(face)

    ! Rotate face until first local index is min global
    DO 
       ! Check if first node matches global min node
       IF (globalMin == globalFace(1)) EXIT
       
       globalFace(1:3) = CSHIFT(globalFace,1)
    END DO
    ! Assign new local numbering
    Tetra % NodeIndexes(face) = globalFace(1:3)

    ! Face 3 now contains global max
    face = getTetraFaceMap(3)
    globalFace(1:3) = Tetra % NodeIndexes(face)
    ! Rotate face until last local index is max global
    DO
       ! Check if last node matches global max node
       IF (globalMax == globalFace(3)) EXIT

       globalFace(1:3) = CSHIFT(globalFace,1)
    END DO
    ! Assign new local numbering
    Tetra % NodeIndexes(face) = globalFace(1:3)

    ! Set AC tetra type
    IF (Tetra % NodeIndexes(2) < Tetra % NodeIndexes(3)) THEN
       Tetra % PDefs % TetraType = 1
    ELSE IF (Tetra % NodeIndexes(3) < Tetra % NodeIndexes(2)) THEN
       Tetra % PDefs % TetraType = 2
    ELSE 
       CALL Fatal('MeshUtils::ConvertToACTetra','Corrupt element type')
    END IF
   
  END SUBROUTINE ConvertToACTetra


!------------------------------------------------------------------------------
!>     Assign local number of edge to given boundary element. Also copies all 
!>     p element attributes from element edge to boundary edge.
!------------------------------------------------------------------------------
  SUBROUTINE AssignLocalNumber( EdgeElement, Element, Mesh,NoPE )
!------------------------------------------------------------------------------
    USE PElementMaps, ONLY : getFaceEdgeMap 
    IMPLICIT NONE

    ! Parameters
    TYPE(Mesh_t) :: Mesh            !< Finite element mesh containing faces and edges.
    TYPE(Element_t), POINTER :: EdgeElement  !< Edge element to which assign local number
    TYPE(Element_t), POINTER :: Element      !< Bulk element with some global numbering to use to assign local number
    LOGICAL, OPTIONAL :: NoPE
!------------------------------------------------------------------------------
    ! Local variables

    INTEGER i,j,n,edgeNumber, numEdges, bMap(4)
    TYPE(Element_t), POINTER :: Edge
    LOGICAL :: EvalPE

    EvalPE = .TRUE.
    IF(PRESENT(NoPE)) EvalPE = .NOT.NoPE

    ! Get number of points, edges or faces
    numEdges = 0
    SELECT CASE (Element % TYPE % DIMENSION)
    CASE (1)
      RETURN
    CASE (2)
       numEdges = Element % TYPE % NumberOfEdges
    CASE (3)   
       numEdges = Element % TYPE % NumberOfFaces
    CASE DEFAULT
       WRITE (*,*) 'MeshUtils::AssignLocalNumber, Unsupported dimension:', Element % TYPE % DIMENSION
       RETURN
    END SELECT

    ! For each edge or face in element try to find local number
    DO edgeNumber=1, numEdges
       ! If edges have not been created, stop search. This should not happen, actually.
       IF (.NOT. ASSOCIATED(Element % EdgeIndexes)) THEN
          ! EdgeElement % localNumber = 0
          RETURN
       END IF

       Edge => GetElementEntity(Element,edgeNumber,Mesh)

       ! Edge element not found. This should not be possible, unless there
       ! is an error in the mesh read in process..
       IF (.NOT. ASSOCIATED(Edge)) THEN
          CALL Warn('MeshUtils::AssignLocalNumber','Edge element not found')
          ! EdgeElement % localNumber = 0
          RETURN
       END IF

       n = 0
       ! For each element node
       DO i=1, Edge % TYPE % NumberOfNodes
          ! For each node in edge element
          DO j=1, EdgeElement % TYPE % NumberOfNodes
             ! If edge and edgeelement node match increment counter
             IF (Edge % NodeIndexes(i) == EdgeElement % NodeIndexes(j)) n = n + 1
          END DO
       END DO

       ! If all nodes are on boundary, edge was found
       IF (n == EdgeElement % TYPE % NumberOfNodes) THEN
          IF(EvalPE) THEN
              EdgeElement % PDefs % localNumber = edgeNumber
              EdgeElement % PDefs % LocalParent => Element
          END IF

          ! Change ordering of global nodes to match that of element
          bMap = getElementBoundaryMap( Element, edgeNumber )
          DO j=1,n
          	EdgeElement % NodeIndexes(j) = Element % NodeIndexes(bMap(j))
	  END DO

          ! Copy attributes of edge element to boundary element
          ! Misc attributes
          IF(EvalPE) THEN
            EdgeElement % PDefs % isEdge = Edge % PDefs % isEdge
          
          ! Gauss points
            EdgeElement % PDefs % GaussPoints = Edge % PDefs % GaussPoints

          ! Element p
            EdgeElement % PDefs % P = Edge % PDefs % P
          END IF
          
          !(and boundary bubble dofs)
          EdgeElement % BDOFs = MAX(EdgeElement % BDOFs, Edge % BDOFs)


          ! If this boundary has edges copy edge indexes
          IF (ASSOCIATED(Edge % EdgeIndexes)) THEN
             ! Allocate element edges to element
             n = Edge % TYPE % NumberOfEdges
             bmap(1:4) = getFaceEdgeMap( Element, edgeNumber )
             
             IF ( ASSOCIATED( EdgeElement % EdgeIndexes) ) THEN
                DEALLOCATE( EdgeElement % EdgeIndexes )
             END IF
             
             CALL AllocateVector( EdgeElement % EdgeIndexes, n )
             ! Copy edges from edge to boundary edge
             DO i=1,n
                EdgeElement % EdgeIndexes(i) = Element % EdgeIndexes(bmap(i))
             !    EdgeElement % EdgeIndexes(i) = Element % EdgeIndexes(i)
             END DO
          END IF
          
          ! Edge fields copied and local edge found so return
          RETURN
       END IF
    END DO

    ! If we are here local number not found
    IF(.NOT.ASSOCIATED(EdgeElement % PDefs % LocalParent)) &
      CALL Warn('MeshUtils::AssignLocalNumber','Unable to find local edge')
    ! EdgeElement % localNumber = 1
  CONTAINS

    FUNCTION GetElementEntity(Element, which, Mesh) RESULT(Entity)
      IMPLICIT NONE

      TYPE(Element_t), POINTER :: Element, Entity 
      INTEGER :: which
      TYPE(Mesh_t) :: Mesh

      NULLIFY(Entity)
      ! Switch by element dimension
      SELECT CASE (Element % TYPE % DIMENSION)
         CASE (2)
            Entity => Mesh % Edges( Element % EdgeIndexes(which))
         CASE (3)
            Entity => Mesh % Faces( Element % FaceIndexes(which))
         CASE DEFAULT
            WRITE (*,*) 'AssignLocalNumber::GetElementEntity: Unsupported dimension'
            RETURN
      END SELECT
    END FUNCTION GetElementEntity
  END SUBROUTINE AssignLocalNumber
    

!------------------------------------------------------------------------------
!>     Based on element degrees of freedom, return the sum of element
!>     degrees of freedom.
!------------------------------------------------------------------------------
  FUNCTION getElementMaxDOFs( Mesh, Element ) RESULT(dofs)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh        !< Finite element mesh
    TYPE(Element_t), POINTER :: Element  !< Element to get maximum dofs for
    INTEGER :: dofs                      !< maximum number of dofs for Element
!------------------------------------------------------------------------------

    TYPE(ELement_t), POINTER :: Edge, Face
    INTEGER :: i, edgeDofs, faceDofs
    
    ! Get sum of edge dofs if any
    edgeDofs = 0
    IF (ASSOCIATED(Element % EdgeIndexes)) THEN
       DO i=1, Element % TYPE % NumberOfEdges
          Edge => Mesh % Edges(Element % EdgeIndexes(i))
          edgeDofs = edgeDofs + Edge % BDOFs
       END DO
    END IF

    ! Get sum of face dofs if any
    faceDofs = 0
    IF (ASSOCIATED(Element % FaceIndexes)) THEN
       DO i=1, Element % TYPE % NumberOfFaces
          Face => Mesh % Faces(Element % FaceIndexes(i))
          faceDofs = faceDofs + Face % BDOFs
       END DO
    END IF

    ! Get sum of all dofs in element
    dofs = Element % TYPE % NumberOfNodes + &
         edgeDofs + faceDofs + Element % BDOFs
  END FUNCTION getElementMaxDOFs




!------------------------------------------------------------------------------
!> Creates a permutation table for bodies or boundaries using a free chosen string
!> as mask. The resulting permutation is optimized in order, if requested. The
!> subroutine is intended to help in saving boundary data in an ordered manner,
!> but it can find other uses as well. Currently the implementation is limited
!> to normal Lagrangian elements.
!------------------------------------------------------------------------------
  SUBROUTINE MakePermUsingMask( Model,Solver,Mesh,MaskName, &
      OptimizeBW, Perm, LocalNodes, MaskOnBulk, RequireLogical, &
      ParallelComm, BreakLoop )
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Mesh_t)   :: Mesh
    TYPE(SOlver_t) :: Solver
    INTEGER :: LocalNodes
    LOGICAL :: OptimizeBW
    INTEGER, POINTER :: Perm(:)
    CHARACTER(LEN=*) :: MaskName
    LOGICAL, OPTIONAL :: MaskOnBulk
    LOGICAL, OPTIONAL :: RequireLogical
    LOGICAL, OPTIONAL :: ParallelComm
    LOGICAL, OPTIONAL :: BreakLoop
!------------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm(:), Neighbours(:)
    INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:)
    TYPE(ListMatrix_t), POINTER :: ListMatrix(:)
    INTEGER :: t,i,j,k,l,m,k1,k2,n,p,q,e1,e2,f1,f2,This,bf_id,nn,t0,ii(ParEnv % PEs)
    INTEGER :: ierr, status(MPI_STATUS_SIZE), NewDofs
    LOGICAL :: Flag, Found, FirstRound, MaskIsLogical, Hit, Parallel
    LOGICAL, ALLOCATABLE :: IsNeighbour(:)
    INTEGER :: Indexes(30), ElemStart, ElemFin, Width, BreakNode
    TYPE(ListMatrixEntry_t), POINTER :: CList, Lptr
    TYPE(Element_t), POINTER :: CurrentElement,Elm
    REAL(KIND=dp) :: MinDist, Dist
!------------------------------------------------------------------------------

    IF(PRESENT(ParallelComm)) THEN
      Parallel = ParallelComm
    ELSE
      Parallel = ParEnv % PEs > 1
    END IF

    ! First check if there are active elements for this mask
    IF( PRESENT( MaskOnBulk ) ) MaskOnBulk = .FALSE.
    IF( PRESENT( RequireLogical ) ) THEN
      MaskIsLogical = RequireLogical
    ELSE
      MaskIsLogical = .FALSE.
    END IF

    IF(.NOT. ASSOCIATED( Perm ) ) THEN
      ALLOCATE( Perm( Mesh % NumberOfNodes ) )
      Perm = 0
    END IF

    ElemStart = HUGE(ElemStart) 
    ElemFin = 0     
    DO l = 1, Model % NumberOfBodyForces
       IF( MaskIsLogical ) THEN
         Hit = ListGetLogical( Model % BodyForces(l) % Values,MaskName,Found) 
       ELSE
         Hit = ListCheckPresent( Model % BodyForces(l) % Values,MaskName)
       END IF 
       IF( Hit ) THEN
          ElemStart = 1
          ElemFin = Mesh % NumberOfBulkElements
          IF( PRESENT( MaskOnBulk ) ) MaskOnBulk = .TRUE.
          EXIT
       END IF
    END DO
    DO l = 1, Model % NumberOfBCs
       IF( MaskIsLogical ) THEN
         Hit = ListGetLogical(Model % BCs(l) % Values,MaskName,Found )
       ELSE
         Hit = ListCheckPresent(Model % BCs(l) % Values,MaskName )
       END IF
       IF( Hit ) THEN
          ElemStart = MIN( ElemStart, Mesh % NumberOfBulkElements + 1)
          ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
          EXIT
       END IF
    END DO

    IF( ElemFin - ElemStart <= 0) THEN
       LocalNodes = 0
       RETURN
    END IF

    k = 0
    Perm = 0
    FirstRound = .TRUE.
    BreakNode = 0
    t0 = 0
    
    ! Loop over the active elements
    ! 1st round initial numbering is given
    ! 2nd round a list matrix giving all the connections is created

100 DO t=ElemStart, ElemFin
       
       CurrentElement => Mesh % Elements(t)
       
       Hit = .FALSE.
       IF(t <= Mesh % NumberOfBulkElements) THEN
          l = CurrentElement % BodyId
	  bf_id = ListGetInteger( Model % Bodies(l) % Values, 'Body Force',Found)
	  IF( bf_id>0 ) THEN
            IF( MaskIsLogical ) THEN
              Hit = ListGetLogical( Model % BodyForces(bf_id) % Values, MaskName, Found )
            ELSE
              Hit = ListCheckPresent( Model % BodyForces(bf_id) % Values, MaskName )
            END IF
	  END IF 
       ELSE
          DO l=1, Model % NumberOfBCs
            IF ( Model % BCs(l) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
            IF( MaskIsLogical ) THEN
              Hit = ListGetLogical(Model % BCs(l) % Values,MaskName, Found ) 
            ELSE
              Hit = ListCheckPresent(Model % BCs(l) % Values,MaskName ) 
            END IF
            EXIT
          END DO
       END IF       
       IF( .NOT. Hit ) CYCLE       
       
       n = CurrentElement % TYPE % NumberOfNodes
       Indexes(1:n) = CurrentElement % NodeIndexes(1:n)
       
       IF( FirstRound ) THEN
         ! Just plainly create the permutation
         DO i=1,n
             j = Indexes(i)
             IF ( Perm(j) == 0 ) THEN
                k = k + 1
                Perm(j) = k
             END IF
          END DO
        ELSE
          ! Create the list matrix for the connectivity in order to minimize the bandwidth
          DO i=1,n
             k1 = Perm(Indexes(i))
             IF ( k1 <= 0 ) CYCLE
             DO j=1,n
                k2 = Perm(Indexes(j))
                IF ( k2 <= 0 ) CYCLE
                IF( k1 == BreakNode .OR. k2 == BreakNode ) THEN
                  IF( t0 == 0 ) t0 = t
                  IF( t0 /= t ) THEN
                    PRINT *,'breaking connection between:',k1,k2
                    CYCLE
                  END IF
                END IF
                Lptr => List_GetMatrixIndex( ListMatrix,k1,k2 )
             END DO
          END DO
       END IF
    END DO
    LocalNodes = k

    ! In parallel case, detect nodes which are shared with another partition
    ! which may not have an element on this boundary
    ! Code borrowed from CommunicateLinearSystemTag
    !------------------------------------------------------------------------------
    IF( Parallel ) THEN

      ALLOCATE( IsNeighbour(ParEnv % PEs), fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )

      nn = MeshNeighbours(Mesh, IsNeighbour)
      nn = 0
      ineigh = 0
      DO i=0, ParEnv % PEs-1
        k = i+1
        IF(i==ParEnv % myPE) CYCLE
        IF(.NOT. IsNeighbour(k) ) CYCLE
        nn = nn + 1
        fneigh(nn) = k
        ineigh(k) = nn
      END DO

      n = COUNT(Perm > 0 .AND. Mesh % ParallelInfo % NodeInterface)
      ALLOCATE( s_e(n, nn ), r_e(n) )

      CALL CheckBuffer( nn*3*n )

      ii = 0
      DO i=1, Mesh % NumberOfNodes
        IF(Perm(i) > 0 .AND. Mesh % ParallelInfo % NodeInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k == ParEnv % MyPE ) CYCLE
            k = k + 1
            k = ineigh(k)
            IF ( k> 0) THEN
              ii(k) = ii(k) + 1
              s_e(ii(k),k) = Mesh % ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
        END IF
      END DO

      DO i=1, nn
        j = fneigh(i)
        CALL MPI_BSEND( ii(i),1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD,ierr )
        IF( ii(i) > 0 ) THEN
          CALL MPI_BSEND( s_e(1:ii(i),i),ii(i),MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,ierr )
        END IF
      END DO

      NewDofs = 0

      DO i=1, nn
        j = fneigh(i)
        CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
        IF ( n>0 ) THEN
          IF( n>SIZE(r_e)) THEN
            DEALLOCATE(r_e)
            ALLOCATE(r_e(n))
          END IF

          CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
          DO j=1,n
            k = SearchNode( Mesh % ParallelInfo, r_e(j), Order=Mesh % ParallelInfo % Gorder )
            IF ( k>0 ) THEN
              IF(.NOT. Perm(k) > 0) THEN
                NewDofs = NewDofs + 1
                Perm(k) = LocalNodes + NewDofs
              END IF
            END IF
          END DO
        END IF
      END DO
      DEALLOCATE(s_e, r_e )

      LocalNodes = LocalNodes + NewDofs
    END IF

    ! Don't optimize bandwidth for parallel cases
    IF( Parallel .OR. .NOT. OptimizeBW ) RETURN

    IF(FirstRound) THEN
       ! Allocate space 
       NULLIFY( ListMatrix )
       ListMatrix => List_AllocateMatrix(LocalNodes)
       FirstRound = .FALSE.

       ! Find the node in the lower left corner at give it the 1st index
       ! since it will probably determine the 1st index
       MinDist = HUGE(MinDist)
       DO i=1,SIZE(Perm)
          IF( Perm(i) <= 0) CYCLE
          Dist = Mesh % Nodes % x(i) + Mesh % Nodes % y(i) + Mesh % Nodes % z(i)
          IF(Dist < MinDist) THEN
             MinDist = Dist
             j = i
          END IF
       END DO

       ! Find the 1st node and swap it with the lower corner
       DO i=1,SIZE(Perm)
          IF( Perm(i) == 1) EXIT
       END DO       
       Perm(i) = Perm(j)
       Perm(j) = 1

       ! Minimizing the bandwidth of a closed loop is impossible.
       ! So let us break the loop on one node. 
       IF(PRESENT(BreakLoop)) THEN
         IF(BreakLoop) BreakNode = 1
       END IF
       
       GOTO 100
    END IF

!------------------------------------------------------------------------------

    ALLOCATE( InvPerm(LocalNodes) )
    InvPerm = 0
    DO i=1,SIZE(Perm)
       IF (Perm(i)>0) InvPerm(Perm(i)) = i
    END DO

    ! The bandwidth optimization for lines results to perfectly ordered 
    ! permutations. If there is only one line the 1st node should be the 
    ! lower left corner.

    Flag = .TRUE.
    Width = OptimizeBandwidth( ListMatrix, Perm, InvPerm, &
        LocalNodes, Flag, Flag, MaskName )

    ! We really only need the permutation, as there will be no matrix equation
    ! associated with it.
    DEALLOCATE( InvPerm )
    CALL List_FreeMatrix( LocalNodes, ListMatrix )

!------------------------------------------------------------------------------
  END SUBROUTINE MakePermUsingMask
!------------------------------------------------------------------------------




!------------------------------------------------------------------------
!> Find a point in the mesh structure
!> There are two strategies:
!> 1) Recursive where the same routine is repeated with sloppier criteria
!> 2) One-sweep strategy where the best hit is registered and used if of 
!>    acceptable accuracy. 
!> There are two different epsilons that control the search. One for the 
!> rough test in absolute coordinates and another one for the more accurate
!> test in local coordinates.   
!-------------------------------------------------------------------------
  FUNCTION PointInMesh(Solver, GlobalCoords, LocalCoords, HitElement, &
      CandElement, ExtInitialize ) RESULT ( Hit )
        
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: GlobalCoords(3), LocalCoords(3)
    TYPE(Element_t), POINTER :: HitElement 
    TYPE(Element_t), POINTER, OPTIONAL :: CandElement
    LOGICAL, OPTIONAL :: ExtInitialize
    LOGICAL :: Hit
!-------------------------------------------------------------------------
    LOGICAL :: Initialize, Allocated = .FALSE., Stat, DummySearch, &
        MaskExists, Found, IsRecursive
    INTEGER :: i,j,k,n,bf_id,dim,mini
    REAL(KIND=dp) :: u,v,w,dist,mindist,MinLocalCoords(3)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Quadrant_t), POINTER, SAVE :: RootQuadrant =>NULL(), LeafQuadrant
    REAL(kind=dp) :: BoundingBox(6), eps2, eps1 = 1d-3, GlobalEps, LocalEps
    CHARACTER(LEN=MAX_NAME_LEN) :: MaskName


    SAVE :: Allocated, ElementNodes, DummySearch, Mesh, MaskName, MaskExists, &
        GlobalEps, LocalEps, IsRecursive


    IF( PRESENT( ExtInitialize ) ) THEN
      Initialize = ExtInitialize
    ELSE
      Initialize = .NOT. Allocated 
    END IF

    IF( Initialize ) THEN
      Mesh => Solver % Mesh
      n = Mesh % MaxElementNodes
      IF( Allocated ) THEN
        DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
      END IF
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n))
      Allocated = .TRUE.

      IsRecursive = ListGetLogical( CurrentModel % Simulation,&
          'Interpolation Search Recursive',Stat )
!      IF(.NOT. Stat ) IsRecursive = .TRUE.

      LocalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Local Epsilon', Stat )
      IF(.NOT. stat) LocalEps = 1.0d-10

      GlobalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Global Epsilon', Stat ) 
      IF(.NOT. stat) THEN
        IF( IsRecursive ) THEN
          GlobalEps = 2.0d-10
        ELSE
          GlobalEps = 1.0d-4
        END IF
      END IF

      DummySearch = ListGetLogical( CurrentModel % Simulation,&
          'Interpolation Search Dummy',Stat )

      MaskName = ListGetString( CurrentModel % Simulation,&
          'Interpolation Search Mask',MaskExists )

      IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
        CALL FreeQuadrantTree( Mesh % RootQuadrant )
        Mesh % RootQuadrant => NULL()
      END IF
    END IF
      

    !-----------------------------------------------
    ! Create the octree search structure, if needed 
    !-----------------------------------------------
    IF ( .NOT. ( DummySearch .OR.  ASSOCIATED( Mesh % RootQuadrant ) ) ) THEN
      BoundingBox(1) = MINVAL( Mesh % Nodes % x )
      BoundingBox(2) = MINVAL( Mesh % Nodes % y )
      BoundingBox(3) = MINVAL( Mesh % Nodes % z )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % x )
      BoundingBox(5) = MAXVAL( Mesh % Nodes % y )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
      
      eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
      BoundingBox(1:3) = BoundingBox(1:3) - eps2
      BoundingBox(4:6) = BoundingBox(4:6) + eps2
      
      CALL BuildQuadrantTree( Mesh,BoundingBox,Mesh % RootQuadrant)
      RootQuadrant => Mesh % RootQuadrant
      IF (.NOT. ASSOCIATED(RootQuadrant) ) THEN
        Hit = .FALSE.
        CALL Warn('PointInMesh','No RootQuadrant associated')
        RETURN
      END IF
    END IF


    Hit = .FALSE.

    ! Check that the previous hit is not hit even now
    !-------------------------------------------------
    IF( PRESENT( CandElement ) ) THEN

      IF( ASSOCIATED(CandElement)) THEN

        CurrentElement => CandElement
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        IF ( PointInElement( CurrentElement, ElementNodes, &
            GlobalCoords, LocalCoords ) ) THEN
          Hit = .TRUE.
          HitElement => CurrentElement
          RETURN
        END IF
      END IF
    END IF


    Eps1 = GlobalEps
    Eps2 = LocalEps


100 IF( DummySearch ) THEN

      mindist = HUGE( mindist ) 
      
      !----------------------------------------------------------
      ! Go through all bulk elements in a dummy search.
      ! This algorithm is mainly here for debugging purposes, or
      ! if just a few nodes need to be searched.
      !----------------------------------------------------------
      DO k=1,Mesh % NumberOfBulkElements
        CurrentElement => Mesh % Elements(k)
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        IF( MaskExists ) THEN
          bf_id = ListGetInteger( CurrentModel % Bodies(CurrentElement % BodyId) % Values, &
              'Body Force', Found )
          IF( .NOT. Found ) CYCLE
          IF(.NOT. ListCheckPresent( CurrentModel % BodyForces(bf_id) % Values,MaskName) ) CYCLE
        END IF

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        Hit = PointInElement( CurrentElement, ElementNodes, &
            GlobalCoords, LocalCoords, Eps1, Eps2, LocalDistance = dist )
        IF( dist < mindist ) THEN
          mini = k
          mindist = dist
        END IF
        IF( Hit ) EXIT
      END DO      
    ELSE
      !-----------------------------------------------
      ! Find the right element using an octree search
      ! This is the preferred algorithms of the two.
      !-----------------------------------------------
      NULLIFY(CurrentElement)
      CALL FindLeafElements(GlobalCoords, Mesh % MeshDim, RootQuadrant, LeafQuadrant)
      IF ( ASSOCIATED(LeafQuadrant) ) THEN
        DO j=1, LeafQuadrant % NElemsInQuadrant
          k = LeafQuadrant % Elements(j)
          CurrentElement => Mesh % Elements(k)
          
          IF( MaskExists ) THEN
            bf_id = ListGetInteger( CurrentModel % Bodies(CurrentElement % BodyId) % Values, &
                'Body Force', Found )
            IF( .NOT. Found ) CYCLE
            IF(.NOT. ListCheckPresent( CurrentModel % BodyForces(bf_id) % Values,MaskName) ) CYCLE
          END IF
          
          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes
                    
          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
          
          Hit = PointInElement( CurrentElement, ElementNodes, &
              GlobalCoords, LocalCoords, Eps1, Eps2, LocalDistance = dist ) 
          IF( dist < mindist ) THEN
            mini = k
            mindist = dist
            MinLocalCoords = LocalCoords
          END IF
          IF( Hit ) EXIT
        END DO
      END IF      
    END IF

    IF( .NOT. Hit ) THEN
      IF( IsRecursive ) THEN
        Eps1 = 10.0 * Eps1
        Eps2 = 10.0 * Eps2
        IF( Eps1 <= 1.0_dp ) GOTO 100
      ELSE
        IF( mindist < Eps1 ) THEN
          CurrentElement => Mesh % Elements(k)
          LocalCoords = MinLocalCoords
          Hit = .TRUE.
        END IF
      END IF
    END IF

    IF( Hit ) HitElement => CurrentElement
    
  END FUNCTION PointInMesh



!--------------------------------------------------------------------------
!> This subroutine finds the structure of an extruded mesh even though it is 
!> given in an unstructured format. The routine may be used by some special
!> solvers that employ the special character of the mesh.
!> The extrusion is found for a given direction and for each node the corresponding 
!> up and down, and thereafter top and bottom node is computed.
!-----------------------------------------------------------------------------
  SUBROUTINE DetectExtrudedStructure( Mesh, Solver, ExtVar, &
      TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer, &
      MidNodePointer, MidLayerExists, NumberOfLayers, NodeLayer, &
      MaskVar )
    
    USE CoordinateSystems
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER, OPTIONAL :: ExtVar
    INTEGER, POINTER, OPTIONAL :: TopNodePointer(:), BotNodePointer(:), &
        UpNodePointer(:), DownNodePointer(:), MidNodePointer(:)
    INTEGER, POINTER, OPTIONAL :: NodeLayer(:)
    INTEGER, OPTIONAL :: NumberOfLayers
    LOGICAL, OPTIONAL :: MidLayerExists
    TYPE(Variable_t), POINTER, OPTIONAL :: MaskVar
!-----------------------------------------------------------------------------
    REAL(KIND=dp) :: Direction(3)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: Tolerance
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    TYPE(Nodes_t), POINTER :: MeshNodes
    INTEGER :: i,j,k,n,ii,jj,dim, nsize, nnodes, elem, TopNodes, BotNodes, Rounds, ActiveDirection, &
	UpHit, DownHit, bc_ind, jmin, jmax
    INTEGER, POINTER :: NodeIndexes(:), MaskPerm(:)
    LOGICAL :: MaskExists, UpActive, DownActive, GotIt, Found, DoCoordTransform
    LOGICAL, POINTER :: TopFlag(:), BotFlag(:)
    REAL(KIND=dp) :: at0, at1, Length, UnitVector(3), Vector(3), Vector2(3), &
        ElemVector(3), DotPro, MaxDotPro, MinDotPro, Eps, MinTop, &
        MaxTop, MinBot, MaxBot
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: TopPointer(:), BotPointer(:), UpPointer(:), DownPointer(:),Layer(:),MidPointer(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CoordTransform
    CHARACTER(*), PARAMETER :: Caller="DetectExtrudedStructure"
   
    CALL Info(Caller,'Determining extruded structure',Level=6)
    at0 = CPUTime()

    DIM = Mesh % MeshDim
    Params => Solver % Values
    
    ActiveDirection = ListGetInteger(Params,'Active Coordinate')
    IF( ActiveDirection < 1 .OR. ActiveDirection > 3 ) THEN
      CALL Fatal('StructuredMeshMapper','Invalid value for Active Coordinate')
    END IF  
    UnitVector = 0.0_dp
    UnitVector(ActiveDirection) = 1.0_dp

    IF( ListGetLogical(Params,'Mapping Original Coordinates',Found ) ) THEN
      MeshNodes => Mesh % NodesOrig
    ELSE
      MeshNodes => Mesh % Nodes
    END IF
    
    IF( ListGetLogical(Params,'Project To Bottom',GotIt) ) &
        UnitVector = -1.0_dp * UnitVector

    WRITE(Message,'(A,3F8.3)') 'Unit vector of direction:',UnitVector
    CALL Info(Caller,Message,Level=8)

    ! Set the dot product tolerance
    !-----------------------------------------------------------------
    Eps = ListGetConstReal( Params,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0d-4

    nnodes = Mesh % NumberOfNodes
    nsize = nnodes

    Var => NULL()
    IF( PRESENT(MaskVar) ) THEN
      Var => MaskVar
    ELSE          
      VarName = ListGetString(Params,'Mapping Mask Variable',GotIt )
      IF(GotIt) THEN
        Var => VariableGet( Mesh % Variables,  VarName )
      END IF
    END IF
    MaskExists = ASSOCIATED(Var)
    IF( MaskExists ) THEN
      ALLOCATE( MaskPerm( SIZE( Var % Perm ) ) )
      MaskPerm = Var % Perm 
      nsize = MAXVAL( MaskPerm ) 
      CALL Info(Caller,'Using variable as mask: '//TRIM(Var % Name),Level=8)
    ELSE
      VarName = ListGetString(Params,'Mapping Mask Name',MaskExists )
      IF( MaskExists ) THEN
        CALL Info(Caller,'Using name as mask: '//TRIM(VarName),Level=8)
        MaskPerm => NULL() 
        CALL MakePermUsingMask( CurrentModel, Solver, Mesh, VarName, &
            .FALSE., MaskPerm, nsize )
        !PRINT *,'nsize:',nsize,SIZE(MaskPerm),MAXVAL(MaskPerm(1:nnodes))
      END IF
    END IF

    IF( MaskExists ) THEN
      CALL Info(Caller,'Applying mask of size: '//TRIM(I2S(nsize)),Level=10)
    ELSE
      CALL Info(Caller,'Applying extrusion on the whole mesh',Level=10)
    END IF 

    CoordTransform = ListGetString(Params,'Mapping Coordinate Transformation',DoCoordTransform )
    IF( DoCoordTransform .OR. MaskExists) THEN
      Var => VariableGet( Mesh % Variables,'Extruded Coordinate')
      IF( ASSOCIATED( Var ) ) THEN
        CALL Info(Caller,'Reusing > Extruded Coordinate < variable',Level=12 )
        Values => Var % Values        
      ELSE
        NULLIFY( Values )
        ALLOCATE( Values( nsize ) )
        Values = 0.0_dp
        IF( MaskExists ) THEN
          CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values, MaskPerm)
        ELSE
          CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values)
        END IF
        Var => VariableGet( Mesh % Variables,'Extruded Coordinate')
      END IF
    ELSE IF( ActiveDirection == 1 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 1')
    ELSE IF( ActiveDirection == 2 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 2')
    ELSE 
      Var => VariableGet( Mesh % Variables,'Coordinate 3')
    END IF	      

    IF( MaskExists .OR. DoCoordTransform) THEN
      DO i=1,Mesh % NumberOfNodes
        j = i
	IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE
        END IF
        Vector(1) = Mesh % Nodes % x(i)
	Vector(2) = Mesh % Nodes % y(i)
	Vector(3) = Mesh % Nodes % z(i)
	IF( DoCoordTransform ) THEN
          CALL CoordinateTransformationNodal( CoordTransform, Vector )
        END IF
        Values(j) = Vector( ActiveDirection )
      END DO
    END IF
    IF( PRESENT( ExtVar ) ) ExtVar => Var
    
    ! Check which direction is active
    !---------------------------------------------------------------------
    UpActive = PRESENT( UpNodePointer) .OR. PRESENT ( TopNodePointer ) 
    DownActive = PRESENT( DownNodePointer) .OR. PRESENT ( BotNodePointer ) 
    
    IF( PRESENT( NumberOfLayers) .OR. PRESENT( NodeLayer ) ) THEN
      UpActive = .TRUE.
      DownActive = .TRUE.
    END IF

    IF(.NOT. (UpActive .OR. DownActive ) ) THEN
      CALL Warn(Caller,'Either up or down direction should be active')
      RETURN
    END IF

    ! Allocate pointers to top and bottom, and temporary pointers up and down
    !------------------------------------------------------------------------
    IF( UpActive ) THEN
      ALLOCATE(TopPointer(nsize),UpPointer(nsize))
      DO i=1,nnodes
        j = i
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE 
        END IF
        TopPointer(j) = i
        UpPointer(j) = i
      END DO
    END IF
    IF( DownActive ) THEN
      ALLOCATE(BotPointer(nsize),DownPointer(nsize))
      DO i=1,nnodes        
        j = i
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE 
        END IF
        BotPointer(j) = i
        DownPointer(j) = i
      END DO
    END IF
    
    CALL Info(Caller,'Determine up and down pointers',Level=15)

    ! Determine the up and down pointers using dot product as criterion
    !-----------------------------------------------------------------
    n = Mesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    
    DO elem = 1,Mesh % NumberOfBulkElements      
      
      Element => Mesh % Elements(elem)
      NodeIndexes => Element % NodeIndexes
      CurrentModel % CurrentElement => Element
      
      n = Element % TYPE % NumberOfNodes
      Nodes % x(1:n) = MeshNodes % x(NodeIndexes)
      Nodes % y(1:n) = MeshNodes % y(NodeIndexes)
      Nodes % z(1:n) = MeshNodes % z(NodeIndexes)
      
      ! This is probably a copy-paste error, I comment it away for time being.   
      ! IF (.NOT. (Element % PartIndex == Parenv % Mype) ) CYCLE

      IF( MaskExists ) THEN
        IF( ANY(MaskPerm(NodeIndexes) == 0) ) CYCLE
      END IF
      
      DO i=1,n
        ii = NodeIndexes(i)
        
        Vector(1) = Nodes % x(i)
	Vector(2) = Nodes % y(i) 
        Vector(3) = Nodes % z(i)
        
 	IF( DoCoordTransform ) THEN
          CALL CoordinateTransformationNodal( CoordTransform, Vector )
        END IF

        MaxDotPro = -1.0_dp
        MinDotPro = 1.0_dp
        
        DO j=i+1,n
          jj = NodeIndexes(j)
          
	  Vector2(1) = Nodes % x(j)
          Vector2(2) = Nodes % y(j)
          Vector2(3) = Nodes % z(j)

	  IF( DoCoordTransform ) THEN
            CALL CoordinateTransformationNodal( CoordTransform, Vector2 )
          END IF
          
          ElemVector = Vector2 - Vector

          Length = SQRT(SUM(ElemVector*ElemVector))
          DotPro = SUM(ElemVector * UnitVector) / Length

          IF( DotPro > MaxDotPro ) THEN
            MaxDotPro = DotPro
            jmax = jj
          END IF
          IF( DotPro < MinDotPro ) THEN
            MinDotPro = DotPro
            jmin = jj
          END IF          
        END DO
          
        IF(MaxDotPro > 1.0_dp - Eps) THEN 
          IF( MaskExists ) THEN
            IF( UpActive ) UpPointer(MaskPerm(ii)) = jmax
            IF( DownActive ) DownPointer(MaskPerm(jmax)) = ii              
          ELSE
            IF( UpActive ) UpPointer(ii) = jmax
            IF( DownActive ) DownPointer(jmax) = ii
          END IF
        END IF
            
        IF(MinDotPro < Eps - 1.0_dp) THEN
          IF( MaskExists ) THEN
            IF( DownActive ) DownPointer(MaskPerm(ii)) = jmin
            IF( UpActive ) UpPointer(MaskPerm(jmin)) = ii
          ELSE
            IF( DownActive ) DownPointer(ii) = jmin
            IF( UpActive ) UpPointer(jmin) = ii              
          END IF
        END IF

      END DO
    END DO
    DEALLOCATE( Nodes % x, Nodes % y,Nodes % z )

    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    CALL Info(Caller,'determine top and bottom pointers',Level=9)

    DO Rounds = 1, nsize
      DownHit = 0
      UpHit = 0
      
      DO i=1,nnodes
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0) CYCLE
          IF( UpActive ) THEN
            j = UpPointer(MaskPerm(i))
            IF( TopPointer(MaskPerm(i)) /= TopPointer(MaskPerm(j)) ) THEN
              UpHit = UpHit + 1
              TopPointer(MaskPerm(i)) = TopPointer(MaskPerm(j))
            END IF
          END IF
          IF( DownActive ) THEN
            j = DownPointer(MaskPerm(i))
            IF( BotPointer(MaskPerm(i)) /= BotPointer(MaskPerm(j)) ) THEN
              DownHit = DownHit + 1
              BotPointer(MaskPerm(i)) = BotPointer(MaskPerm(j))
            END IF
          END IF
        ELSE
          IF( UpActive ) THEN
            j = UpPointer(i)
            IF( TopPointer(i) /= TopPointer(j) ) THEN
              UpHit = UpHit + 1
              TopPointer(i) = TopPointer( j )
            END IF
          END IF
          IF( DownActive ) THEN
            j = DownPointer(i)
            IF( BotPointer(i) /= BotPointer( j ) ) THEN
              DownHit = DownHit + 1
              BotPointer(i) = BotPointer( j )
            END IF
          END IF
        END IF
      END DO
      
      IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
    END DO

    ! The last round is always a check
    Rounds = Rounds - 1
    
    CALL Info(Caller,'Layered structure detected in '//TRIM(I2S(Rounds))//' cycles',Level=9)
    IF( Rounds == 0 ) THEN
      CALL Info(Caller,'Try to increase value for > Dot Product Tolerance < ')
      CALL Fatal(Caller,'Zero rounds implies unsuccessful operation')
    END IF

    ! Compute the number of layers. The Rounds above may in some cases
    ! be too small. Here just one layer is used to determine the number
    ! of layers to save some time.
    !------------------------------------------------------------------
    IF( PRESENT( NumberOfLayers ) ) THEN
      CALL Info(Caller,'Compute number of layers',Level=15)    
      DO i=1,nsize
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0 ) CYCLE
        END IF
        EXIT
      END DO

      j = BotPointer(1)      
      CALL Info(Caller,'Starting from node: '//TRIM(I2S(j)),Level=15)

      NumberOfLayers = 0
      DO WHILE(.TRUE.)
        jj = j 
        IF( MaskExists ) THEN
          jj = MaskPerm(j)
        END IF
        k = UpPointer(jj)
        IF( k == j ) THEN
          EXIT
        ELSE
          NumberOfLayers = NumberOfLayers + 1
          j = k
        END IF
      END DO

      IF( NumberOfLayers < Rounds ) THEN
        WRITE( Message,'(A,I0,A,I0)') 'There seems to be varying number of layers: ',&
            NumberOfLayers,' vs. ',Rounds
        CALL Warn(Caller, Message )
        NumberOfLayers = Rounds
      END IF
      CALL Info(Caller,&
          'Extruded structure layers: '//TRIM(I2S(NumberOfLayers)),Level=6)
    END IF

    
    ! Create layer index if requested
    !------------------------------------------------------------------
    IF( PRESENT( NodeLayer ) ) THEN
      CALL Info(Caller,'creating layer index',Level=9)        

      NULLIFY(Layer)
      ALLOCATE( Layer(nsize) )
      Layer = 1
      IF( MaskExists ) THEN
        WHERE( MaskPerm == 0 ) Layer = 0
        
        DO i=1,nnodes
          IF( MaskPerm(i) == 0 ) CYCLE
          Rounds = 1
          j = BotPointer(MaskPerm(i))
          Layer(MaskPerm(j)) = Rounds
          DO WHILE(.TRUE.)
            k = UpPointer(MaskPerm(j))
            IF( k == j ) EXIT          
            Rounds = Rounds + 1
            j = k
            Layer(MaskPerm(j)) = Rounds
          END DO
        END DO
      ELSE        
        DO i=1,nsize
          Rounds = 1
          j = BotPointer(i)
          Layer(j) = Rounds
          DO WHILE(.TRUE.)
            k = UpPointer(j)
            IF( k == j ) EXIT          
            Rounds = Rounds + 1
            j = k
            Layer(j) = Rounds
          END DO
        END DO
      END IF
        
      NodeLayer => Layer
      WRITE(Message,'(A,I0,A,I0,A)') 'Layer range: [',MINVAL(Layer),',',MAXVAL(Layer),']'
      CALL Info(Caller,Message,Level=6)
      NULLIFY(Layer)
    END IF

    
    IF( PRESENT( MidNodePointer ) ) THEN
      ALLOCATE( MidPointer( nsize ) )
      MidPointer = 0 
      MidLayerExists = .FALSE.

      DO elem = Mesh % NumberOfBulkElements + 1, &       
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements  
        
        Element => Mesh % Elements(elem)
        NodeIndexes => Element % NodeIndexes
        
        DO bc_ind = 1, CurrentModel % NumberOfBCs 
          IF( Element % BoundaryInfo % Constraint == &
              CurrentModel % BCs(bc_ind) % Tag ) THEN
            IF( ListCheckPresent( CurrentModel % BCs(bc_ind) % Values,'Mid Surface') ) THEN
              MidPointer( NodeIndexes ) = NodeIndexes
              MidLayerExists = .TRUE.
            END IF
            EXIT
          END IF
        END DO
      END DO

      IF( MidLayerExists ) THEN
        CALL Info(Caller,'determine mid pointers',Level=15)       
                
        DO Rounds = 1, nsize
          DownHit = 0
          UpHit = 0
          DO i=1,nsize
            IF( MaskExists ) THEN
              IF( MaskPerm(i) == 0) CYCLE
            END IF

            ! We can only start from existing mid pointer
            IF( MidPointer(i) == 0 ) CYCLE
            IF( UpActive ) THEN
              j = UpPointer(i)
              IF( MaskExists ) THEN
                IF( MidPointer(MaskPerm(j)) == 0 ) THEN
                  UpHit = UpHit + 1
                  MidPointer(MaskPerm(j)) = MidPointer(MaskPerm(i))
                END IF
              ELSE
                IF( MidPointer(j) == 0 ) THEN
                  UpHit = UpHit + 1
                  MidPointer(j) = MidPointer(i)
                END IF
              END IF
            END IF
            IF( DownActive ) THEN
              j = DownPointer(i)
              IF( MaskExists ) THEN
                IF( MidPointer(MaskPerm(j)) == 0 ) THEN
                  DownHit = DownHit + 1
                  MidPointer(MaskPerm(j)) = MidPointer(MaskPerm(i))
                END IF           
              ELSE
                IF( MidPointer(j) == 0 ) THEN
                  DownHit = DownHit + 1
                  MidPointer(j) = MidPointer(i)
                END IF
              END IF
            END IF
          END DO
          IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
        END DO

        CALL Info(Caller,&
            'Mid layer structure detected in '//TRIM(I2S(Rounds-1))//' cycles',Level=9)
        MidNodePointer => MidPointer
      ELSE
        DEALLOCATE( MidPointer ) 
        MidNodePointer => NULL()
      END IF
    END IF

  
    ! Count the number of top and bottom nodes, for information only
    !---------------------------------------------------------------
    CALL Info(Caller,'Counting top and bottom nodes',Level=15)        
    IF( UpActive ) THEN
      TopNodes = 0
      MinTop = HUGE( MinTop ) 
      MaxTop = -HUGE( MaxTop )
      DO i=1,nnodes
        IF( MaskExists ) THEN
          j = MaskPerm(i) 
          IF( j == 0 ) CYCLE
          IF(TopPointer(j) == i) THEN
            MinTop = MIN( MinTop, Var % Values(j) )
            MaxTop = MAX( MaxTop, Var % Values(j) )
            TopNodes = TopNodes + 1
          END IF
        ELSE
          IF(TopPointer(i) == i) THEN
            MinTop = MIN( MinTop, Var % Values(i) )
            MaxTop = MAX( MaxTop, Var % Values(i) )
            TopNodes = TopNodes + 1
          END IF
        END IF
      END DO
    END IF

    IF( DownActive ) THEN
      BotNodes = 0
      MinBot = HUGE( MinBot ) 
      MaxBot = -HUGE( MaxBot )
      DO i=1,nnodes
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE
          IF( BotPointer(j) == i) THEN
            MinBot = MIN( MinBot, Var % Values(j))
            MaxBot = MAX( MaxBot, Var % Values(j))
            BotNodes = BotNodes + 1
          END IF
        ELSE          
          IF(BotPointer(i) == i) THEN
            MinBot = MIN( MinBot, Var % Values(i))
            MaxBot = MAX( MaxBot, Var % Values(i))
            BotNodes = BotNodes + 1
          END IF
        END IF
      END DO
    END IF


    ! Return the requested pointer structures, otherwise deallocate
    !---------------------------------------------------------------
    CALL Info(Caller,'Setting pointer structures',Level=15)        
    IF( UpActive ) THEN
      IF( PRESENT( TopNodePointer ) ) THEN
        TopNodePointer => TopPointer 
        NULLIFY( TopPointer )
      ELSE
        DEALLOCATE( TopPointer )
      END IF
      IF( PRESENT( UpNodePointer ) ) THEN
        UpNodePointer => UpPointer 
        NULLIFY( UpPointer )
      ELSE
        DEALLOCATE( UpPointer )
      END IF
    END IF
    IF( DownActive ) THEN
      IF( PRESENT( BotNodePointer ) ) THEN
        BotNodePointer => BotPointer 
        NULLIFY( BotPointer ) 
      ELSE
        DEALLOCATE( BotPointer )
      END IF
      IF( PRESENT( DownNodePointer ) ) THEN
        DownNodePointer => DownPointer 
        NULLIFY( DownPointer ) 
      ELSE
        DEALLOCATE( DownPointer )
      END IF
    END IF

    !---------------------------------------------------------------
    at1 = CPUTime()  
    WRITE(Message,* ) 'Top and bottom pointer init time: ',at1-at0
    CALL Info(Caller,Message,Level=6)
    CALL Info(Caller,&
        'Top and bottom pointer init rounds: '//TRIM(I2S(Rounds)),Level=5)
    IF( UpActive ) THEN
      CALL Info(Caller,'Number of nodes at the top: '//TRIM(I2S(TopNodes)),Level=6)
    END IF
    IF( DownActive ) THEN
      CALL Info(Caller,'Number of nodes at the bottom: '//TRIM(I2S(BotNodes)),Level=6)
    END IF

    IF(DownActive .AND. UpActive ) THEN
      IF(TopNodes /= BotNodes ) THEN
        CALL Fatal(Caller, 'Something wrong: top and bottom node counts differ!')
      END IF
    END IF
    

  CONTAINS
    
    
    !---------------------------------------------------------------
    SUBROUTINE CoordinateTransformationNodal( CoordTransform, R )
      CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform
      REAL(KIND=dp) :: R(3)
      !---------------------------------------------------------------
      REAL(KIND=dp) :: Rtmp(3)
      REAL(KIND=dp), SAVE :: Coeff 
      LOGICAL, SAVE :: Visited = .FALSE.
      

      IF( .NOT. Visited ) THEN
        IF( ListGetLogical( Params,'Angles in Degrees') ) THEN
          Coeff = 180.0_dp / PI
        ELSE
          Coeff = 1.0_dp
        END IF
        Visited = .TRUE.
      END IF
      
      SELECT CASE ( CoordTransform )
        
      CASE('cartesian to cylindrical')
        Rtmp(1) = SQRT( R(1)**2 + R(2)**2)
        Rtmp(2) = Coeff * ATAN2( R(2), R(1)  ) 
        Rtmp(3) = R(3) 
        
      CASE('cylindrical to cartesian')
        Rtmp(1) = COS( R(2) / Coeff ) * R(1)
        Rtmp(2) = SIN( R(2) / Coeff ) * R(1)
        Rtmp(3) = R(3)
        
      CASE DEFAULT
        CALL Fatal('CoordinateTransformationNodal','Unknown transformation: '//TRIM(CoordTransform) )
        
      END SELECT
      
      R = Rtmp

    END SUBROUTINE CoordinateTransformationNodal
   

  END SUBROUTINE DetectExtrudedStructure
 !---------------------------------------------------------------



!--------------------------------------------------------------------------
!> This subroutine finds the structure of an extruded mesh for elements.
!> Otherwise very similar as the DetectExtrudedStructure for nodes.
!> Mesh faces may need to be created in order to determine the up and down
!> pointers.
!-----------------------------------------------------------------------------
  SUBROUTINE DetectExtrudedElements( Mesh, Solver, ExtVar, &
      TopElemPointer, BotElemPointer, UpElemPointer, DownElemPointer, &
      NumberOfLayers, ElemLayer )
    
    USE CoordinateSystems
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER, OPTIONAL :: ExtVar
    INTEGER, POINTER, OPTIONAL :: TopElemPointer(:), BotElemPointer(:), &
        UpElemPointer(:), DownElemPointer(:)
    INTEGER, POINTER, OPTIONAL :: ElemLayer(:)
    INTEGER, OPTIONAL :: NumberOfLayers
!-----------------------------------------------------------------------------
    REAL(KIND=dp) :: Direction(3)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: Tolerance
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(Nodes_t) :: Nodes
    TYPE(Nodes_t), POINTER :: MeshNodes
    INTEGER :: i,j,k,n,ii,jj,dim, nsize, elem, TopNodes, BotNodes, Rounds, ActiveDirection, &
	UpHit, DownHit, bc_ind
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: UpActive, DownActive, GotIt, Found
    LOGICAL, POINTER :: TopFlag(:), BotFlag(:)
    REAL(KIND=dp) :: at0, at1
    REAL(KIND=dp) :: FaceCenter(3),FaceDx(3),Height(2),Eps, MinTop, MaxTop, MinBot, MaxBot, Diam
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: TopPointer(:), BotPointer(:), UpPointer(:), DownPointer(:),Layer(:),MidPointer(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName
    INTEGER :: TestCounter(3),ElementIndex(2)
    CHARACTER(*),PARAMETER :: Caller="DetectExtrudedElements"
         
    CALL Info(Caller,'Determining extruded element structure',Level=6)
    at0 = CPUTime()

    DIM = Mesh % MeshDim

    IF( DIM /= 3 ) THEN
      CALL Fatal(Caller,'Only implemented for 3D cases: '//TRIM(I2S(dim)))
    END IF

    IF( .NOT. ASSOCIATED( Mesh % Faces ) ) THEN
      CALL FindMeshFaces3D( Mesh )
    END IF

    
    Params => Solver % Values
    TestCounter = 0
    
    ActiveDirection = ListGetInteger(Params,'Active Coordinate')
    IF( ActiveDirection < 1 .OR. ActiveDirection > 3 ) THEN
      CALL Fatal(Caller,'Invalid value for Active Coordinate')
    END IF

    IF( ListGetLogical(Params,'Mapping Original Coordinates',Found ) ) THEN
      MeshNodes => Mesh % NodesOrig
    ELSE
      MeshNodes => Mesh % Nodes
    END IF
    
    ! Set the dot product tolerance
    !-----------------------------------------------------------------
    Eps = ListGetConstReal( Params,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0d-1

    nsize = Mesh % NumberOfBulkElements
    CALL Info(Caller,'Detecting extrusion in the mesh using coordinate: '&
        //TRIM(I2S(ActiveDirection)),Level=8)

    IF( ActiveDirection == 1 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 1')
    ELSE IF( ActiveDirection == 2 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 2')
    ELSE 
      Var => VariableGet( Mesh % Variables,'Coordinate 3')
    END IF	      

    IF( PRESENT( ExtVar ) ) ExtVar => Var

    ! Check which direction is active
    !---------------------------------------------------------------------
    UpActive = PRESENT( UpElemPointer) .OR. PRESENT ( TopElemPointer ) 
    DownActive = PRESENT( DownElemPointer) .OR. PRESENT ( BotElemPointer ) 

    IF( PRESENT( NumberOfLayers) .OR. PRESENT( ElemLayer ) ) THEN
      UpActive = .TRUE.
      DownActive = .TRUE.
    END IF

    IF(.NOT. (UpActive .OR. DownActive ) ) THEN
      CALL Warn(Caller,'Either up or down direction should be active')
      RETURN
    END IF

    ! Allocate pointers to top and bottom, and temporary pointers up and down
    !------------------------------------------------------------------------
    IF( UpActive ) THEN
      ALLOCATE(TopPointer(nsize),UpPointer(nsize))
      DO i=1,nsize
        TopPointer(i) = i
        UpPointer(i) = i
      END DO
    END IF
    IF( DownActive ) THEN
      ALLOCATE(BotPointer(nsize),DownPointer(nsize))
      DO i=1,nsize
        BotPointer(i) = i
        DownPointer(i) = i
      END DO
    END IF

    CALL Info(Caller,'determine up and down pointers',Level=15)

    ! Determine the up and down pointers using dot product as criterion
    !-----------------------------------------------------------------
    n = Mesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    
    DO elem = 1,Mesh % NumberOfFaces 

      Element => Mesh % Faces(elem)
      NodeIndexes => Element % NodeIndexes
      CurrentModel % CurrentElement => Element

      n = Element % TYPE % NumberOfNodes
      Nodes % x(1:n) = MeshNodes % x(NodeIndexes)
      Nodes % y(1:n) = MeshNodes % y(NodeIndexes)
      Nodes % z(1:n) = MeshNodes % z(NodeIndexes)

      IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right ) ) CYCLE
      
      FaceCenter(1) = SUM( Nodes % x(1:n) ) / n
      FaceCenter(2) = SUM( Nodes % y(1:n) ) / n
      FaceCenter(3) = SUM( Nodes % z(1:n) ) / n

      FaceDx(1) = SUM( ABS( Nodes % x(1:n) - FaceCenter(1) ) ) 
      FaceDx(2) = SUM( ABS( Nodes % y(1:n) - FaceCenter(2) ) ) 
      FaceDx(3) = SUM( ABS( Nodes % z(1:n) - FaceCenter(3) ) ) 
      
      Diam = SQRT( SUM( FaceDx**2 ) )

      ! This is not a face that separates extruded elements
      IF( FaceDx(ActiveDirection) > Eps * Diam ) CYCLE      

      TestCounter(1) = TestCounter(1) + 1      
      
      DO k = 1, 2
        IF( k == 1 ) THEN
          Parent => Element % BoundaryInfo % Left
        ELSE
          Parent => Element % BoundaryInfo % Right
        END IF
        IF( .NOT. ASSOCIATED( Parent ) ) CYCLE
               
        n = Parent % TYPE % NumberOfNodes
        NodeIndexes => Parent % NodeIndexes        

        ElementIndex(k) = Parent % ElementIndex
        Height(k) = SUM( Var % Values(NodeIndexes) ) / n
      END DO      

      IF( Height(1) > Height(2) ) THEN
        IF( UpActive ) UpPointer(ElementIndex(2)) = ElementIndex(1)
        IF( DownActive ) DownPointer(ElementIndex(1)) = ElementIndex(2)
      ELSE
        IF( UpActive ) UpPointer(ElementIndex(1)) = ElementIndex(2)
        IF( DownActive ) DownPointer(ElementIndex(2)) = ElementIndex(1)
      END IF
    END DO  
        
    DEALLOCATE( Nodes % x, Nodes % y,Nodes % z )

    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    CALL Info(Caller,'determine top and bottom pointers',Level=9)

    DO Rounds = 1, nsize
      DownHit = 0
      UpHit = 0
      DO i=1,nsize
        IF( UpActive ) THEN
          j = UpPointer(i)
          IF( TopPointer(i) /= TopPointer( j ) ) THEN
            UpHit = UpHit + 1
            TopPointer(i) = TopPointer( j )
          END IF
        END IF
        IF( DownActive ) THEN
          j = DownPointer(i)
          IF( BotPointer(i) /= BotPointer( j ) ) THEN
	    DownHit = DownHit + 1
            BotPointer(i) = BotPointer( j )
          END IF
        END IF
      END DO
      CALL Info(Caller,'Hits in determining structure: '//TRIM(I2S(UpHit+DownHit)),Level=10)
      IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
    END DO
    ! The last round is always a check
    Rounds = Rounds - 1


    WRITE( Message,'(A,I0,A)') 'Layered elements detected in ',Rounds,' cycles'
    CALL Info(Caller,Message,Level=9)
    IF( Rounds == 0 ) THEN
      CALL Info(Caller,'Try to increase value for > Dot Product Tolerance < ')
      CALL Fatal(Caller,'Zero rounds implies unsuccessful operation')
    END IF


    ! Compute the number of layers. The Rounds above may in some cases 
    ! be too small. Here just one layer is used to determine the number
    ! of layers to save some time.
    !------------------------------------------------------------------
    IF( PRESENT( NumberOfLayers ) ) THEN
      CALL Info(Caller,'Compute number of layers',Level=15)    

      ! We start from any bottom row entry
      j = BotPointer(1)
      
      NumberOfLayers = 0
      DO WHILE(.TRUE.)
        k = UpPointer(j)

        IF( k == j ) THEN
          EXIT
        ELSE
          NumberOfLayers = NumberOfLayers + 1
          j = k
        END IF
      END DO      

      IF( NumberOfLayers < Rounds ) THEN
        WRITE( Message,'(A,I0,A,I0)') 'There seems to be varying number of layers: ',&
            NumberOfLayers,' vs. ',Rounds
        CALL Warn(Caller, Message )
        NumberOfLayers = Rounds
      END IF
      CALL Info(Caller,'Extruded structure layers: '//TRIM(I2S(NumberOfLayers)),Level=6)
    END IF

    
    ! Create layer index if requested
    !------------------------------------------------------------------
    IF( PRESENT( ElemLayer ) ) THEN
      CALL Info(Caller,'creating layer index',Level=9)        

      NULLIFY(Layer)
      ALLOCATE( Layer(nsize) )
      Layer = 1
      
      DO i=1,nsize
        Rounds = 1
        j = BotPointer(i)
        Layer(j) = Rounds
        DO WHILE(.TRUE.)
          k = UpPointer(j)
          IF( k == j ) EXIT          
          Rounds = Rounds + 1
          j = k
          Layer(j) = Rounds
        END DO
      END DO
      
      ElemLayer => Layer
      WRITE(Message,'(A,I0,A,I0,A)') 'Layer range: [',MINVAL(Layer),',',MAXVAL(Layer),']'
      CALL Info(Caller,Message,Level=6)
      NULLIFY(Layer)
    END IF

  
    ! Count the number of top and bottom elements, for information only
    !---------------------------------------------------------------
    CALL Info(Caller,'Counting top and bottom elements',Level=15)        
    IF( UpActive ) THEN
      TopNodes = 0
      MinTop = HUGE( MinTop ) 
      MaxTop = -HUGE( MaxTop )
      DO i=1,nsize
        IF(TopPointer(i) == i) THEN
          MinTop = MIN( MinTop, Var % Values(i) )
          MaxTop = MAX( MaxTop, Var % Values(i) )
          TopNodes = TopNodes + 1
        END IF
      END DO
      CALL Info(Caller,'Number of top elements: '//TRIM(I2S(TopNodes)),Level=9)
    END IF

    IF( DownActive ) THEN
      BotNodes = 0
      MinBot = HUGE( MinBot ) 
      MaxBot = -HUGE( MaxBot )
      DO i=1,nsize
        IF(BotPointer(i) == i) THEN
          MinBot = MIN( MinBot, Var % Values(i))
          MaxBot = MAX( MaxBot, Var % Values(i))
          BotNodes = BotNodes + 1
        END IF
      END DO
    END IF


    ! Return the requested pointer structures, otherwise deallocate
    !---------------------------------------------------------------
    CALL Info(Caller,'Setting pointer structures',Level=15)        
    IF( UpActive ) THEN
      IF( PRESENT( TopElemPointer ) ) THEN
        TopElemPointer => TopPointer 
        NULLIFY( TopPointer )
      ELSE
        DEALLOCATE( TopPointer )
      END IF
      IF( PRESENT( UpElemPointer ) ) THEN
        UpElemPointer => UpPointer 
        NULLIFY( UpPointer )
      ELSE
        DEALLOCATE( UpPointer )
      END IF
    END IF
    IF( DownActive ) THEN
      IF( PRESENT( BotElemPointer ) ) THEN
        BotElemPointer => BotPointer 
        NULLIFY( BotPointer ) 
      ELSE
        DEALLOCATE( BotPointer )
      END IF
      IF( PRESENT( DownElemPointer ) ) THEN
        DownElemPointer => DownPointer 
        NULLIFY( DownPointer ) 
      ELSE
        DEALLOCATE( DownPointer )
      END IF
    END IF

    !---------------------------------------------------------------
    at1 = CPUTime()  
    WRITE(Message,'(A,ES12.3)') 'Top and bottom pointer init time: ',at1-at0
    CALL Info(Caller,Message,Level=6)

    CALL Info(Caller,'Top and bottom pointer init rounds: '//TRIM(I2S(Rounds)),Level=8)

    IF( UpActive ) THEN
      CALL Info(Caller,'Number of elements at the top: '//TRIM(I2S(TopNodes)),Level=8)
    END IF
    IF( DownActive ) THEN
      CALL Info(Caller,'Number of elements at the bottom: '//TRIM(I2S(BotNodes)),Level=8)
    END IF

    IF(DownActive .AND. UpActive ) THEN
      IF(TopNodes /= BotNodes ) THEN
        CALL Fatal(Caller, 'Something wrong: top and bottom element counts differ!')
      END IF
    END IF
    

  END SUBROUTINE DetectExtrudedElements
 !---------------------------------------------------------------


  SUBROUTINE StoreOriginalCoordinates(Mesh)
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER CONTIG :: NewCoords(:)
    INTEGER :: n

    IF( ASSOCIATED( Mesh % NodesOrig ) ) THEN
      CALL Info('StoreOriginalCoordinates','Original coordinates already stored')
    END IF

    n = SIZE( Mesh % Nodes % x )    
    NULLIFY( NewCoords )
    ALLOCATE( NewCoords(3*n) )

    ALLOCATE( Mesh % NodesOrig ) 
    Mesh % NodesOrig % x => NewCoords(1:n)
    Mesh % NodesOrig % y => NewCoords(n+1:2*n)
    Mesh % NodesOrig % z => NewCoords(2*n+1:3*n)

    Mesh % NodesOrig % x = Mesh % Nodes % x
    Mesh % NodesOrig % y = Mesh % Nodes % y
    Mesh % NodesOrig % z = Mesh % Nodes % z

    Mesh % NodesMapped => Mesh % Nodes

    CALL Info('StoreOriginalCoordinates','Original coordinates stored',Level=6)
    
  END SUBROUTINE StoreOriginalCoordinates

    
   
  !----------------------------------------------------------------
  !> Maps coordinates from the original nodes into a new coordinate
  !> system while optionally maintaining the original coordinates. 
  !> Note that this may be called 
  !---------------------------------------------------------------
  SUBROUTINE CoordinateTransformation( Mesh, CoordTransform, Params, &
      IrreversibleTransformation )
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL, OPTIONAL :: IrreversibleTransformation
    !---------------------------------------------------------------   
    REAL(KIND=dp) :: R0(3),R1(3),Coeff,Rad0
    LOGICAL :: Irreversible,FirstTime,Reuse,UpdateNodes,Found
    REAL(KIND=dp), POINTER :: x0(:),y0(:),z0(:),x1(:),y1(:),z1(:)
    REAL(KIND=dp), POINTER CONTIG :: NewCoords(:)
    INTEGER :: i,j,k,n,Mode
    TYPE(Variable_t), POINTER :: Var

    ! The coordinate transformation may either be global for all the solvers
    ! and this overrides the original nodes permanently. 
    ! Or it can be a solver specific transformation which saves the initial 
    ! coordinates. 
    CALL Info('CoordinateTransformation','Starting')

    IF(.NOT. ASSOCIATED(Mesh) ) THEN
      CALL Fatal('CoordinateTransformation','Mesh not associated!')
    END IF

    IF( PRESENT( IrreversibleTransformation ) ) THEN
      Irreversible = IrreversibleTransformation
    ELSE
      Irreversible = .FALSE.
    END IF

    n = Mesh % NumberOfNodes 

    x0 => Mesh % Nodes % x
    y0 => Mesh % Nodes % y
    z0 => Mesh % Nodes % z
    
    IF( Irreversible ) THEN
      UpdateNodes = .TRUE.
      ! Map to the same nodes
      x1 => Mesh % Nodes % x
      y1 => Mesh % Nodes % y
      z1 => Mesh % Nodes % z
    ELSE
      ReUse = ListGetLogical(Params,'Coordinate Transformation Reuse',Found ) 
      FirstTime = .NOT. ASSOCIATED( Mesh % NodesMapped )
      IF( FirstTime ) THEN
        ALLOCATE( Mesh % NodesMapped )
        NULLIFY( NewCoords )
        ALLOCATE( NewCoords(3*n) )
        NewCoords = 0.0_dp
        Mesh % NodesMapped % x => NewCoords(1:n)
        Mesh % NodesMapped % y => NewCoords(n+1:2*n)
        Mesh % NodesMapped % z => NewCoords(2*n+1:3*n)
        ! Mesh % NodesMapped % x => NewCoords(1::3)
        ! Mesh % NodesMapped % y => NewCoords(2::3)
        ! Mesh % NodesMapped % z => NewCoords(3::3)
      ELSE
        IF( n /= SIZE(Mesh % NodesMapped % x) ) THEN
          CALL Fatal('CoordinateTransformation','Sizes of original and mapped mesh differ!')
        END IF
      END IF

      IF( CoordTransform == 'previous' ) THEN
        IF( FirstTime ) THEN
          CALL Fatal('CoordinateTransformation','One cannot reuse unexisting transformation!')
        END IF
        ReUse = .TRUE.
      END IF

      ! Note that if many solvers reutilize the same coordinates then they must 
      ! also have the same coordinate mapping. 
      !------------------------------------------------------------------------
      UpdateNodes = FirstTime .OR. .NOT. ReUse 
      ! Map different nodes if the original ones are kept
      x1 => Mesh % NodesMapped % x
      y1 => Mesh % NodesMapped % y
      z1 => Mesh % NodesMapped % z      

      IF( FirstTime ) THEN
        IF( ListGetLogical(Params,'Coordinate Transformation Save',Found ) ) THEN
          CALL Info('CoordinateTranformation',&
              'Creating variables for > Transformed Coordinate < ')
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 1',1,x1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 2',1,y1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 3',1,z1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate',3,NewCoords)
        END IF
      END IF
    END IF
      
    IF( UpdateNodes ) THEN
      IF( ListGetLogical( Params,'Coordinate Transformation Use Degrees',Found) ) THEN
        Coeff = 180.0_dp / PI
        CALL Info('CoordinateTranformation','Using degrees for angles')
      ELSE
        Coeff = 1.0_dp
      END IF

      Rad0 = ListGetConstReal( Params,'Coordinate Transformation Radius',Found )
  
      SELECT CASE ( CoordTransform ) 
        
      CASE('cartesian to polar')
        Mode = 1
      CASE('cartesian to cylindrical')
        Mode = 1
      CASE('polar to cartesian')
        Mode = -1
      CASE('cylindrical to cartesian')
        Mode = -1
        
      CASE DEFAULT
        CALL Fatal('CoordinateTransformation','Unknown transformation: '//TRIM(CoordTransform) )
        
      END SELECT

      DO i=1,n    
        R0(1) = x0(i)
        R0(2) = y0(i)
        R0(3) = z0(i)
        
        IF( Mode == 1 ) THEN
          R1(1) = Rad0 + SQRT( R0(1)**2 + R0(2)**2)
          R1(2) = Coeff * ATAN2( R0(2), R0(1)  ) 
          R1(3) = R0(3)    
       
        ELSE IF( Mode == -1 ) THEN
          R1(1) = COS( R0(2) / Coeff ) * ( R0(1) + Rad0 )
          R1(2) = SIN( R0(2) / Coeff ) * ( R0(1) + Rad0 )
          R1(3) = R0(3)          
        END IF

        x1(i) = R1(1)
        y1(i) = R1(2)
        z1(i) = R1(3)

      END DO
    END IF

    IF( .NOT. Irreversible ) THEN
      Mesh % NodesOrig => Mesh % Nodes
      Mesh % Nodes => Mesh % NodesMapped

      Var => VariableGet( CurrentModel % Variables,'Coordinate 1')
      Var % Values => Mesh % Nodes % x

      Var => VariableGet( CurrentModel % Variables,'Coordinate 2')
      Var % Values => Mesh % Nodes % y

      Var => VariableGet( CurrentModel % Variables,'Coordinate 3')
      Var % Values => Mesh % Nodes % z
    END IF

    CALL Info('CoordinateTransformation','All done',Level=8)

  END SUBROUTINE CoordinateTransformation
!---------------------------------------------------------------

  

!---------------------------------------------------------------
!> Return back to the original coordinate system. 
!---------------------------------------------------------------
  SUBROUTINE BackCoordinateTransformation( Mesh, DeleteTemporalMesh )
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: DeleteTemporalMesh
!---------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var

    IF( PRESENT( DeleteTemporalMesh ) ) THEN
      IF( DeleteTemporalMesh ) THEN
        DEALLOCATE( Mesh % NodesMapped % x, &
            Mesh % NodesMapped % y, &
            Mesh % NodesMapped % z ) 
        DEALLOCATE( Mesh % NodesMapped )
      END IF
    END IF

    IF( .NOT. ASSOCIATED( Mesh % NodesOrig ) ) THEN
      CALL Fatal('BackCoordinateTransformation','NodesOrig not associated')
    END IF

    Mesh % Nodes => Mesh % NodesOrig

    Var => VariableGet( CurrentModel % Variables,'Coordinate 1')
    Var % Values => Mesh % Nodes % x
    
    Var => VariableGet( CurrentModel % Variables,'Coordinate 2')
    Var % Values => Mesh % Nodes % y

    Var => VariableGet( CurrentModel % Variables,'Coordinate 3')
    Var % Values => Mesh % Nodes % z

  END SUBROUTINE BackCoordinateTransformation
!---------------------------------------------------------------


 
  !> Find the node closest to the given coordinate. 
  !> The linear search only makes sense for a small number of points. 
  !> Users include saving routines of pointwise information. 
  !-----------------------------------------------------------------
  FUNCTION ClosestNodeInMesh(Mesh,Coord,MinDist) RESULT ( NodeIndx )
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coord(3)
    REAL(KIND=dp), OPTIONAL :: MinDist
    INTEGER :: NodeIndx

    REAL(KIND=dp) :: Dist2,MinDist2,NodeCoord(3)
    INTEGER :: i

    MinDist2 = HUGE( MinDist2 ) 

    DO i=1,Mesh % NumberOfNodes
      
      NodeCoord(1) = Mesh % Nodes % x(i)
      NodeCoord(2) = Mesh % Nodes % y(i)
      NodeCoord(3) = Mesh % Nodes % z(i)
    
      Dist2 = SUM( ( Coord - NodeCoord )**2 )
      IF( Dist2 < MinDist2 ) THEN
        MinDist2 = Dist2
        NodeIndx = i  
      END IF
    END DO
    
    IF( PRESENT( MinDist ) ) MinDist = SQRT( MinDist2 ) 

  END FUNCTION ClosestNodeInMesh


  !> Find the element that owns or is closest to the given coordinate. 
  !> The linear search only makes sense for a small number of points. 
  !> Users include saving routines of pointwise information. 
  !-------------------------------------------------------------------
  FUNCTION ClosestElementInMesh(Mesh, Coords) RESULT ( ElemIndx )

    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coords(3)
    INTEGER :: ElemIndx

    REAL(KIND=dp) :: Dist,MinDist,LocalCoords(3)
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: k,l,n,istat
    REAL(KIND=dp) :: ParallelHits,ParallelCands
    LOGICAL :: Hit

    n = Mesh % MaxElementNodes
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), STAT=istat)
    IF( istat /= 0 ) CALL Fatal('ClosestElementInMesh','Memory allocation error') 	
    ElemIndx = 0
    MinDist = HUGE( MinDist ) 
    Hit = .FALSE.
    l = 0
    
    ! Go through all bulk elements and look for hit in each element.
    ! Linear search makes only sense for a small number of nodes
    DO k=1,Mesh % NumberOfBulkElements

      Element => Mesh % Elements(k)
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      
      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
      
      Hit = PointInElement( Element, ElementNodes, &
          Coords, LocalCoords, LocalDistance = Dist )
      IF( Dist < MinDist ) THEN
        MinDist = Dist
        l = k
      END IF
      IF( Hit ) EXIT
    END DO
    
    ! Count the number of parallel hits
    !-----------------------------------------------------------------------
    IF( Hit ) THEN
      ParallelHits = 1.0_dp
    ELSE
      ParallelHits = 0.0_dp
    END IF
    ParallelHits = ParallelReduction( ParallelHits )
    
    ! If there was no proper hit go through the best candidates so far and 
    ! see if they would give a acceptable hit
    !----------------------------------------------------------------------
    IF( ParallelHits < 0.5_dp ) THEN	  

      ! Compute the number of parallel candidates
      !------------------------------------------
      IF( l > 0 ) THEN
        ParallelCands = 1.0_dp
      ELSE
        ParallelCands = 0.0_dp
      END IF
      ParallelCands = ParallelReduction( ParallelCands ) 

      IF( l > 0 ) THEN
        Element => Mesh % Elements(l)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

        ! If there are more than two competing parallel hits then use more stringent conditions
        ! since afterwards there is no way of deciding which one was closer.
        !--------------------------------------------------------------------------------------
        IF( ParallelCands > 1.5_dp ) THEN
          Hit = PointInElement( Element, ElementNodes, &
              Coords, LocalCoords, GlobalEps = 1.0d-3, LocalEps=1.0d-4 )	
        ELSE
          Hit = PointInElement( Element, ElementNodes, &
              Coords, LocalCoords, GlobalEps = 1.0_dp, LocalEps=0.1_dp )	
        END IF
      END IF
    END IF

    IF( Hit ) ElemIndx = l

    IF( ParallelHits < 0.5_dp ) THEN
      IF( Hit ) THEN
        ParallelHits = 1.0_dp
      ELSE
        ParallelHits = 0.0_dp
      END IF
      ParallelHits = ParallelReduction( ParallelHits )
      IF( ParallelHits < 0.5_dp ) THEN
        WRITE( Message, * ) 'Coordinate not found in any of the elements!',Coords
        CALL Warn( 'ClosestElementInMesh', Message )
      END IF
    END IF

    DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
 
  END FUNCTION ClosestElementInMesh



!---------------------------------------------------------------
!> This find two fixing nodes for each coordinate direction
!> The indexes are returned in order: x1 x2 y1 y2 z1 z2.
!---------------------------------------------------------------
  SUBROUTINE FindRigidBodyFixingNodes(Solver,FixingDofs,MaskPerm)
!------------------------------------------------------------------------------
    USE GeneralUtils

    TYPE(Solver_t) :: Solver
    INTEGER, OPTIONAL :: FixingDofs(0:)
    INTEGER, OPTIONAL :: MaskPerm(:)

!---------------------------------------------------------------

    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: MaskExists,FixBestDirection,FoundBetter, GotIt
    INTEGER :: i,j,k,l,ind,n,dim,dir,nsize,Sweep,MaxSweep,DirBest
    INTEGER :: PosMeasureIndex, NegMeasureIndex, FixingNodes(0:6)
    LOGICAL, ALLOCATABLE :: ForbiddenNodes(:)
    REAL(KIND=dp), POINTER :: Parray(:,:)
    REAL(KIND=dp) :: Normal(3), Tangent1(3), Tangent2(3), Coord(3), &
        SumCoord(3), AveCoord(3), Weights(3), RefScore, Score, &
        PosMeasure, NegMeasure, OffLineCoeff, DirDistance, &
        InLine, OffLine, Dist, MinDist, InLineMeasure, ScoreLimit
    CHARACTER(LEN=MAX_NAME_LEN) :: Method
!---------------------------------------------------------------

    CALL Info('FindRigidBodyFixingNodes','Starting',Level=6)

    Mesh => Solver % Mesh
    dim = Mesh % MeshDim 
    
    ALLOCATE( ForbiddenNodes(Mesh % NumberOfNodes) )
    CALL DetermineForbiddenNodes( )
    nsize = COUNT(.NOT. ForbiddenNodes) 

!   PRINT *,'Number of allowed Nodes:',nsize

    ! Find the center from the average of node positions
    !-----------------------------------------------------------
    SumCoord = 0.0_dp
    DO i=1,Mesh % NumberOfNodes
      IF( ForbiddenNodes( i ) ) CYCLE
      
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
    
      SumCoord = SumCoord + Coord
    END DO
    AveCoord = SumCoord / nsize


    ! Find the node closest to center and make that the new center
    !--------------------------------------------------------------
    MinDist = HUGE( MinDist ) 

    DO i=1,Mesh % NumberOfNodes
      IF( ForbiddenNodes( i ) ) CYCLE
      
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
    
      Dist = SUM( ( Coord - AveCoord )**2 )
      IF( Dist < MinDist ) THEN
        MinDist = Dist
        k = i  
      END IF
    END DO

    AveCoord(1) = Mesh % Nodes % x(k)
    AveCoord(2) = Mesh % Nodes % y(k)
    AveCoord(3) = Mesh % Nodes % z(k)
    IF(PRESENT(FixingDOFs)) FixingDOFs(0)=k
    

!   PRINT *,'AveCoord:',AveCoord

    ! Parameters of the search
    !-----------------------------------------------------------

    OffLineCoeff = ListGetConstReal( Solver % Values,'Fixing Nodes Off Line Coefficient',GotIt)
    IF(.NOT. GotIt) OffLineCoeff = 1.0_dp

    ScoreLimit = ListGetConstReal( Solver % Values,'Fixing Nodes Limit Score',GotIt)
    IF(.NOT. GotIt) ScoreLimit = 0.99_dp

    FixBestDirection = ListGetLogical( Solver % Values,'Fixing Nodes Axis Freeze',GotIt)

    Parray => ListGetConstRealArray( Solver % Values,'Fixing Nodes Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal = 0.0_dp
      Normal(1) = 1.0
    END IF
    Normal = Normal / SQRT( SUM( Normal ** 2) )      
    CALL TangentDirections( Normal,Tangent1,Tangent2 )
    
    ! Find the fixing nodes by looping over all nodes
    !-----------------------------------------------------------
    DirDistance = 0.0_dp
    DirBest = 0
    DO dir = 1, dim
      
      ! Use the three principal directions as the weight
      !-------------------------------------------------
      IF( dir == 1 ) THEN
        Weights = Normal
      ELSE IF( dir == 2 ) THEN
        Weights = Tangent1
      ELSE 
        Weights = Tangent2
      END IF
      
      PosMeasure = 0.0_dp
      PosMeasureIndex = 0
      NegMeasure = 0.0_dp
      NegMeasureIndex = 0


      ! Choose the nodes within the cones in the given three directions
      !---------------------------------------------------------------
      DO i=1,Mesh % NumberOfNodes
        IF( ForbiddenNodes( i ) ) CYCLE
        
        Coord(1) = Mesh % Nodes % x(i) 
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)
        
        Coord = Coord - AveCoord
        Dist = SQRT( SUM( Coord ** 2 ) )
 
        ! Signed distance in in-line direction
        InLine = SUM( Coord * Weights )
        
        ! Distance in off-line direction 
        OffLine = SQRT( Dist**2 - InLine**2 )
        
        ! This defines a cone within which nodes are accepted
        InLineMeasure = ABS( InLine ) - OffLineCoeff * OffLine 
        IF( InLineMeasure < 0.0_dp ) CYCLE
        
        IF( InLine < 0.0_dp ) THEN
          IF( InLineMeasure > NegMeasure ) THEN
            NegMeasure = InLineMeasure
            NegMeasureIndex = i
          END IF
        ELSE           
          IF( InLineMeasure > PosMeasure ) THEN
            PosMeasure = InLineMeasure 
            PosMeasureIndex = i
          END IF
        END IF      
      END DO
      
      FixingNodes(2*dir-1) = NegMeasureIndex
      FixingNodes(2*dir) = PosMeasureIndex      

      IF( NegMeasureIndex > 0 .AND. PosMeasureIndex > 0 ) THEN
        IF( PosMeasure + NegMeasure > DirDistance ) THEN
          DirDistance = PosMeasure + NegMeasure
          DirBest = dir
        END IF
      END IF

    END DO


 
    ! To be on the safe side check that no node is used twice
    ! However, do not break the best direction
    !-----------------------------------------------------------------------------------
    DO i=1,2*dim
      DO j=1,2*dim
        IF( FixBestDirection ) THEN
          IF( j == 2*DirBest-1 .OR. j == 2*DirBest ) CYCLE
        END IF        
        IF( FixingNodes(j) == FixingNodes(i) ) FixingNodes(j) = 0
      END DO
    END DO


    ! Go through the fixing nodes one-by-one and set the node so that the harmonic sum
    ! is minimized. This means that small distances are hopefully eliminated. 
    !-----------------------------------------------------------------------------------
    MaxSweep = ListGetInteger( Solver % Values,'Fixing Nodes Search Loops',GotIt)
    DO Sweep = 0,MaxSweep
      FoundBetter = .FALSE.
      DO j=1,2*dim 
        RefScore = FixingNodesScore(j,FixingNodes(j)) 

        ! The first round set the unfixed nodes
        IF( Sweep == 0 ) THEN
!         PRINT *,'Initial Score:',j,RefScore
          IF( FixingNodes(j) /= 0 ) CYCLE
        END IF

        ! Fir the best direction because otherwise there are too 
        ! many moving parts.
        IF( FixBestDirection ) THEN
          IF( j == 2*DirBest-1 .OR. j == 2*DirBest ) CYCLE
        END IF

        RefScore = FixingNodesScore(j,FixingNodes(j)) 

        DO i=1,Mesh % NumberOfNodes
          IF( ForbiddenNodes(i) ) CYCLE
          Score = FixingNodesScore(j,i)
          IF( Score < ScoreLimit * RefScore ) THEN
            RefScore = Score 
            FixingNodes(j) = i            
            FoundBetter = .TRUE.
          END IF
        END DO
      END DO
      IF(.NOT. FoundBetter ) EXIT
    END DO

    DO j=1,2*dim
      RefScore = FixingNodesScore(j,FixingNodes(j)) 
!     PRINT *,'Final Score:',j,RefScore
    END DO

    ! Output the selected nodes
    !-----------------------------------------------------------------------------------
    DO i=1,2*dim
      j = FixingNodes(i)
      WRITE(Message,'(A,I0,3ES10.2)') 'Fixing Node: ',j,&
          Mesh % Nodes % x( j ), &
          Mesh % Nodes % y( j ), &
          Mesh % Nodes % z( j ) 
      CALL Info('FindRigidBodyFixingNodes',Message,Level=6)
      IF( PRESENT( FixingDofs ) ) FixingDofs(i) = j     
    END DO

    DEALLOCATE( ForbiddenNodes )


  CONTAINS

    !> Find the nodes that are either on interface, boundary or do not belong to the field.
    !-----------------------------------------------------------------------------------
    SUBROUTINE DetermineForbiddenNodes()

      TYPE(Element_t), POINTER :: Element
      LOGICAL, POINTER :: ig(:)
      INTEGER :: t
      
      ! Mark all interface nodes as forbidden nodes
      !-----------------------------------------------
      IF( ParEnv % PEs > 1 ) THEN
        ig => Mesh % ParallelInfo % NodeInterface
        ForbiddenNodes = ig(1:Mesh % NumberOfNodes)
      END IF

      ! Mark all nodes on boundary elements as forbidden nodes
      !--------------------------------------------------------
      DO t=Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements

        Element => Mesh % Elements( t )
        ForbiddenNodes( Element % NodeIndexes ) = .TRUE.
      END DO

      ! If mask exists then add all nodes not in mask to forbidden nodes
      !-----------------------------------------------------------------
      IF( PRESENT( MaskPerm) ) THEN
        DO i=1,Mesh % NumberOfNodes
          IF( MaskPerm(i) == 0 ) ForbiddenNodes(i) = .TRUE.
        END DO
      END IF
      
    END SUBROUTINE DetermineForbiddenNodes


    !> Give a value of goodness to the chosen fixing node.
    !-----------------------------------------------------------------------------------
    FUNCTION FixingNodesScore(direction,cand) RESULT ( Score )

      INTEGER :: direction, cand
      INTEGER :: i,j
      REAL(KIND=dp) :: Score

      REAL(KIND=dp) :: x0(3), x1(3), Dist

      IF( cand == 0 ) THEN
        Score = HUGE( Score ) 
        RETURN
      END IF

      Score = 0.0_dp
      x0(1) = Mesh % Nodes % x( cand )
      x0(2) = Mesh % Nodes % y( cand )
      x0(3) = Mesh % Nodes % z( cand )

      DO i=1,2*dim
        IF( i == direction ) CYCLE
        j = FixingNodes( i )

        ! Do not measure distance to unset nodes!
        IF( j == 0 ) CYCLE

        ! This would lead to division by zero later on
        IF( cand == j ) THEN
          Score = HUGE( Score ) 
          RETURN
        END IF

        x1(1) = Mesh % Nodes % x( j )
        x1(2) = Mesh % Nodes % y( j )
        x1(3) = Mesh % Nodes % z( j )

        Dist = SQRT( SUM( (x0 - x1 ) ** 2 ) )
        Score = Score + 1 / Dist
      END DO

    END FUNCTION FixingNodesScore


!------------------------------------------------------------------------------
  END SUBROUTINE FindRigidBodyFixingNodes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Create a 1D mesh, may be used in 1D outlet conditions, for example.
!------------------------------------------------------------------------------
  FUNCTION CreateLineMesh( Params ) RESULT( Mesh )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params 
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    INTEGER :: i, j, k, n, NoNodes, NoElements, ActiveDirection, Order, BodyId, ne
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t),POINTER :: elmt
    REAL(KIND=dp) :: MeshVector(3), Length, Coord(3)
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshName
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    
!------------------------------------------------------------------------------
    Mesh => NULL()
    IF ( .NOT. ASSOCIATED( Params ) ) RETURN
    Mesh => AllocateMesh()

    CALL Info('CreateLineMesh','Creating 1D mesh on-the-fly')

!   Read in the parameters defining a uniform 1D mesh
!--------------------------------------------------------------    
    Order = ListGetInteger( Params,'1D Element Order',Found,minv=1,maxv=2)
    NoElements = ListGetInteger( Params,'1D Number Of Elements',minv=1)
    Length = ListGetConstReal( Params,'1D Mesh Length',Found)
    IF(.NOT. Found) Length = 1.0_dp
    ActiveDirection = ListGetInteger( Params,'1D Active Direction',Found,minv=-3,maxv=3)
    IF(.NOT.Found) ActiveDirection = 1
    BodyId = ListGetInteger( Params,'1D Body Id',Found,minv=1)
    IF(.NOT. Found) BodyId = 1
    MeshName = ListGetString( Params,'1D Mesh Name',Found)
    IF(.NOT. Found) MeshName = '1d_mesh'
    
    Mesh % Name = MeshName
    Mesh % OutputActive = .FALSE.

!   Compute the resulting mesh parameters
!--------------------------------------------------------------
    ne = Order + 1
    NoNodes = NoElements + 1 + NoElements * (Order - 1)    
    MeshVector = 0.0_dp
    MeshVector( ABS( ActiveDirection ) ) = 1.0_dp
    IF( ActiveDirection < 0 ) MeshVector = -MeshVector
    MeshVector = MeshVector * Length
    
!   Define nodal coordinates
!   -------------------------------
    CALL AllocateVector( Mesh % Nodes % x, NoNodes )
    CALL AllocateVector( Mesh % Nodes % y, NoNodes )
    CALL AllocateVector( Mesh % Nodes % z, NoNodes )

    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z

    ALLOCATE( w(0:NoNodes-1) )
    
    CALL UnitSegmentDivision( w, NoNodes-1, Params )
    
    DO i=1, NoNodes
      Coord = MeshVector * w(i-1)

      x(i) = Coord(1)
      y(i) = Coord(2)
      z(i) = Coord(3)
    END DO
    

!   Define elements
!   -------------------------------
    CALL AllocateVector( Mesh % Elements, NoElements )

    Elmt => GetElementType( 200 + ne )

    DO i=1,NoElements
      Element => Mesh % Elements(i)      
      Element % TYPE => Elmt
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()     
      Element % ElementIndex = i

      CALL AllocateVector( Element % NodeIndexes, ne )
      Element % Ndofs = ne ! TO DO: This is not consistent for "Element = n:N", with N>1

      Element % NodeIndexes(1) = (i-1)*Order + 1
      Element % NodeIndexes(2) = i*Order + 1

      DO j=3,ne
        Element % NodeIndexes(j) = (i-1)*Order + j-1
      END DO
      
      Element % BodyId = BodyId
      Element % PartIndex = ParEnv % myPE
    END DO
    
!   Update new mesh node count:
!   ---------------------------

    Mesh % NumberOfNodes = NoNodes
    Mesh % Nodes % NumberOfNodes = NoNodes
    Mesh % NumberOfBulkElements = NoElements
    Mesh % MaxElementNodes = ne
    Mesh % MaxElementDOFs = ne
    Mesh % MeshDim = 1

    CALL SetMeshMaxDOFs(Mesh)

    
    WRITE(Message,'(A,I0)') 'Number of elements created: ',NoElements
    CALL Info('CreateLineMesh',Message)

    WRITE(Message,'(A,I0)') 'Number of nodes created: ',NoNodes
    CALL Info('CreateLineMesh',Message)
 
    CALL Info('CreateLineMesh','All done')

  END FUNCTION CreateLineMesh

  !Creates a regular 2D mesh of 404 elements
  !The resulting mesh has no boundary elements etc for now
  !Should only be used for e.g. mesh to mesh interpolation
  FUNCTION CreateRectangularMesh(Params) RESULT(Mesh)

!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    REAL(KIND=dp) :: min_x, max_x, min_y, max_y, dx, dy
    INTEGER :: i, j, k, n, counter, nnx, nny, nex, ney, &
         NoNodes, NoElements, col, row
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t),POINTER :: elmt
    REAL(KIND=dp) :: MeshVector(3), Length, Coord(3)
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshName, FuncName="CreateRectangularMesh"

!------------------------------------------------------------------------------
    Mesh => NULL()
    IF ( .NOT. ASSOCIATED( Params ) ) RETURN
    Mesh => AllocateMesh()

    CALL Info(FuncName,'Creating 2D mesh on-the-fly')

    !Get parameters from valuelist
    min_x = ListGetConstReal(Params, "Grid Mesh Min X",UnfoundFatal=.TRUE.)
    max_x = ListGetConstReal(Params, "Grid Mesh Max X",UnfoundFatal=.TRUE.)
    min_y = ListGetConstReal(Params, "Grid Mesh Min Y",UnfoundFatal=.TRUE.)
    max_y = ListGetConstReal(Params, "Grid Mesh Max Y",UnfoundFatal=.TRUE.)
    dx    = ListGetConstReal(Params, "Grid Mesh dx",UnfoundFatal=.TRUE.)
    dy    = ListGetConstReal(Params, "Grid Mesh dy",Found)
    IF(.NOT. Found) dy = dx

    IF(max_x <= min_x .OR. max_y <= min_y .OR. dx <= 0.0_dp .OR. dy <= 0.0_dp) &
         CALL Fatal(FuncName, "Bad Grid Mesh parameters!")

    !number of nodes in x and y direction (and total)
    nnx = FLOOR((max_x - min_x) / dx) + 1
    nny = FLOOR((max_y - min_y) / dy) + 1
    NoNodes = nnx * nny

    !number of elements in x and y direction (and total)
    nex = nnx - 1
    ney = nny - 1
    NoElements = nex * ney


!   Define nodal coordinates
!   -------------------------------
    CALL AllocateVector( Mesh % Nodes % x, NoNodes )
    CALL AllocateVector( Mesh % Nodes % y, NoNodes )
    CALL AllocateVector( Mesh % Nodes % z, NoNodes )
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z

    z = 0.0_dp !2D

    !Define node positions
    counter = 0
    DO i=1,nnx
      DO j=1,nny
        counter = counter + 1
        x(counter) = min_x + (i-1)*dx
        y(counter) = min_y + (j-1)*dy
      END DO
    END DO

!   Define elements
!   -------------------------------
    CALL AllocateVector( Mesh % Elements, NoElements )

    Elmt => GetElementType( 404 )

    DO i=1,NoElements
      Element => Mesh % Elements(i)
      Element % TYPE => Elmt
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()
      Element % ElementIndex = i
      CALL AllocateVector( Element % NodeIndexes, 4 )
      Element % Ndofs = 4 ! TO DO: This is not consistent for "Element = n:N", with N>1

      col = MOD(i-1,ney)
      row = (i-1)/ney

      !THIS HERE NEEDS FIXED!!!!!
      Element % NodeIndexes(1) = (row * nny) + col + 1
      Element % NodeIndexes(2) = (row * nny) + col + 2
      Element % NodeIndexes(4) = ((row+1) * nny) + col + 1
      Element % NodeIndexes(3) = ((row+1) * nny) + col + 2

      Element % BodyId = 1
      Element % PartIndex = ParEnv % myPE
    END DO

!   Update new mesh node count:
!   ---------------------------

    Mesh % NumberOfNodes = NoNodes
    Mesh % Nodes % NumberOfNodes = NoNodes
    Mesh % NumberOfBulkElements = NoElements
    Mesh % MaxElementNodes = 4
    Mesh % MaxElementDOFs = 4
    Mesh % MeshDim = 2

  END FUNCTION CreateRectangularMesh

  SUBROUTINE ElmerMeshToDualGraph(Mesh, DualGraph, UseBoundaryMesh)
    IMPLICIT NONE

    TYPE(Mesh_t) :: Mesh
    TYPE(Graph_t) :: DualGraph
    LOGICAL, OPTIONAL :: UseBoundaryMesh

    TYPE(Element_t), POINTER :: Element, Elements(:)

    ! MESH DATA
    ! Mesh (CRS format)
    INTEGER, ALLOCATABLE :: eptr(:), eind(:)
    INTEGER :: nelem
    ! Vertex to element map (CRS format)
    INTEGER, ALLOCATABLE :: vptr(:), vind(:)
    INTEGER :: nvertex

    ! WORK ARRAYS
    ! Pointers to vertex-element maps of the current element
    INTEGER, ALLOCATABLE :: ptrli(:), ptrti(:)
    ! Neighbour indices
    INTEGER, ALLOCATABLE :: neighind(:)
    ! ARRAY MERGE: map for merge
    INTEGER, ALLOCATABLE :: wrkmap(:)

    TYPE :: IntTuple_t
      INTEGER :: i1, i2
    END type IntTuple_t

    TYPE(IntTuple_t), ALLOCATABLE :: wrkheap(:)

    ! OpenMP thread block leads for work division
    INTEGER, ALLOCATABLE :: thrblk(:)
    ! Work indices
    INTEGER, ALLOCATABLE :: wrkind(:), wrkindresize(:)
    INTEGER :: nwrkind

    ! Variables
    INTEGER :: i, dnnz, eid, nl, nli, nti, nn, nv, nthr, &
            te, thrli, thrti, vli, vti, TID, allocstat
    INTEGER :: mapSizePad, maxNodesPad, neighSizePad
    LOGICAL :: Boundary

    INTEGER, PARAMETER :: HEAPALG_THRESHOLD = 24

    CALL Info('ElmerMeshToDualGraph','Creating a dual graph for the mesh',Level=8)

    Boundary = .FALSE.
    IF (Present(UseBoundaryMesh)) Boundary = UseBoundaryMesh

    ! Pointers to mesh data
    IF (.NOT. Boundary) THEN
       nelem = Mesh % NumberOfBulkElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements
    ELSE
       nelem = Mesh % NumberOfBoundaryElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements(&
            Mesh % NumberOfBulkElements+1:Mesh % NumberOfBulkElements+nelem)
    END IF

    ! Initialize dual mesh size and number of nonzeroes
    DualGraph % n = nelem
    dnnz = 0

    ! Copy mesh to CRS structure
    ALLOCATE(eptr(nelem+1), eind(nelem*Mesh % MaxElementNodes), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate mesh structure!')

    eptr(1)=1 ! Fortran numbering
    DO i=1, nelem
      Element => Elements(i)
      nl = Element % TYPE % NumberOfNodes
      nli = eptr(i) ! Fortran numbering
      nti = nli+nl-1
      eind(nli:nti) = Element % NodeIndexes(1:nl) ! Fortran numbering
      eptr(i+1) = nli+nl
    END DO

    ! Construct vertex to element list (in serial!)
    CALL VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)

    ! Allocate pointers to dual mesh
    ALLOCATE(DualGraph % ptr(nelem+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')

    ! Divide work by number of rows in the vertex graph
    nthr = 1 
    !$ nthr = omp_get_max_threads()

    ! Load balance the actual work done by threads (slow)
    ! CALL ThreadLoadBalanceElementNeighbour(nthr, nelem, eptr, eind, vptr, thrblk)
    CALL ThreadStaticWorkShare(nthr, nelem, thrblk)

    !$OMP PARALLEL SHARED(nelem, nvertex, eptr, eind, &
    !$OMP                 vptr, vind, Mesh, DualGraph, &
    !$OMP                 nthr, thrblk, dnnz) &
    !$OMP PRIVATE(i, eid, nli, nti, nn, nv, vli, vti, te, &
    !$OMP         maxNodesPad, neighSizePad, ptrli, ptrti, &
    !$OMP         wrkheap, wrkmap, neighind, &
    !$OMP         wrkind, nwrkind, wrkindresize, allocstat, &
    !$OMP         mapSizePad, thrli, thrti, TID) NUM_THREADS(nthr) &
    !$OMP DEFAULT(NONE)

    TID = 1
    !$ TID = OMP_GET_THREAD_NUM()+1

    ! Ensure that the vertex to element lists are sorted
    !$OMP DO 
    DO i=1,nvertex
      vli = vptr(i)
      vti = vptr(i+1)-1

      CALL Sort(vti-vli+1, vind(vli:vti))
    END DO
    !$OMP END DO NOWAIT

    ! Allocate work array (local to each thread)
    maxNodesPad = IntegerNBytePad(Mesh % MaxElementNodes, 8)
    neighSizePad = IntegerNBytePad(Mesh % MaxElementNodes*20, 8)

    ! Pointers to vertex maps
    ALLOCATE(neighind(neighSizePad), &
            ptrli(maxNodesPad), ptrti(maxNodesPad), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')
    ! Initialize neighbour indices
    neighind = 0

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      ! With multiple threads, use heap based merge
      ALLOCATE(wrkheap(maxNodesPad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
    ELSE
      ! With a small number of threads, use map -based merge
      mapSizePad = IntegerNBytePad(nelem, 8)
      ALLOCATE(wrkmap(mapSizePad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
      ! Initialize local map
      wrkmap=0
    END IF

    ! Allocate local list for results
    nwrkind = 0
    ALLOCATE(wrkind(nelem/nthr*20), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')

    ! Ensure that all the threads have finished sorting the vertex indices
    !$OMP BARRIER

    ! Get thread indices
    thrli = thrblk(TID)
    thrti = thrblk(TID+1)

    ! For each element
    DO eid=thrli,thrti-1
      nli = eptr(eid)
      nti = eptr(eid+1)-1
      nv = nti-nli+1

      ! Get pointers to vertices related to the nodes of the element
      te = 0
      DO i=nli,nti
        ptrli(i-nli+1)=vptr(eind(i))
        ptrti(i-nli+1)=vptr(eind(i)+1) ! NOTE: This is to make comparison cheaper
        te = te + ptrti(i-nli+1)-ptrli(i-nli+1)
      END DO

      ! Allocate neighind large enough
      IF (SIZE(neighind)<te) THEN
        DEALLOCATE(neighind)
        neighSizePad = IntegerNBytePad(te,8)
        ALLOCATE(neighind(neighSizePad), STAT=allocstat)
        neighind = 0
      END IF

      ! Merge vertex lists (multi-way merge of ordered lists)
      IF (nthr >= HEAPALG_THRESHOLD) THEN
        CALL kWayMergeHeap(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkheap)
      ELSE
        CALL kWayMergeArray(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkmap)
      END IF

      ! Add merged list to final list of vertices
      IF (nn+nwrkind>SIZE(wrkind)) THEN
        ALLOCATE(wrkindresize(MAX(nn+nwrkind,2*SIZE(wrkind))), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                'Unable to allocate local workspace!')
        wrkindresize(1:nwrkind)=wrkind(1:nwrkind)
        DEALLOCATE(wrkind)
        CALL MOVE_ALLOC(wrkindresize, wrkind)
      END IF
      wrkind(nwrkind+1:nwrkind+nn) = neighind(1:nn)
      nwrkind = nwrkind + nn

      ! Store number of row nonzeroes
      DualGraph % ptr(eid)=nn
    END DO

    ! Get the global size of the dual mesh
    !$OMP DO REDUCTION(+:dnnz)
    DO i=1,nthr
      dnnz = nwrkind
    END DO
    !$OMP END DO

    ! Allocate memory for dual mesh indices
    !$OMP SINGLE
    ALLOCATE(DualGraph % ind(dnnz), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')
    ! ptr stores row counts, build crs pointers from them
    CALL ComputeCRSIndexes(nelem, DualGraph % ptr)
    !$OMP END SINGLE

    DualGraph % ind(&
            DualGraph % ptr(thrli):DualGraph % ptr(thrti)-1)=wrkind(1:nwrkind)

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      DEALLOCATE(wrkheap, STAT=allocstat)
    ELSE
      DEALLOCATE(wrkmap, STAT=allocstat)
    END IF
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to deallocate local workspace!')
    DEALLOCATE(neighind, ptrli, ptrti, wrkind)

    !$OMP END PARALLEL

    ! Deallocate the rest of memory
    DEALLOCATE(eind, eptr, vptr, vind, thrblk)

    CALL Info('ElmerMeshToDualGraph','Dual graph created with size '//TRIM(I2S(dnnz)),Level=8)


  CONTAINS

    SUBROUTINE VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nelem, nvertex
      INTEGER :: eptr(:), eind(:)
      INTEGER, ALLOCATABLE :: vptr(:), vind(:)

      INTEGER :: i, j, v, eli, eti, ind, tmpi, tmpip, allocstat

      ! Initialize vertex structure (enough storage for nvertex vertices
      ! having eptr(nelem+1) elements)
      ALLOCATE(vptr(nvertex+1), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
              'Vertex allocation failed!')
      vptr = 0

      ! For each element

      ! Compute number of elements attached to each vertex (size of lists)
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        DO j=eli, eti
          vptr(eind(j))=vptr(eind(j))+1
        END DO
      END DO

      ! Compute in-place cumulative sum (row pointers!)
      CALL ComputeCRSIndexes(nvertex, vptr)

      ! Allocate vertex to element lists
      ALLOCATE(vind(vptr(nvertex+1)), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
              'Vertex allocation failed!')

      ! Construct element lists for each vertex
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        ! For each vertex in element
        DO j=eli, eti
          ! Add connection to vertex eind(j)
          ind = eind(j)
          vind(vptr(ind))=i
          vptr(ind)=vptr(ind)+1
        END DO
      END DO

      ! Correct row pointers
      DO i=nvertex,2,-1
        vptr(i)=vptr(i-1)
      END DO
      vptr(1)=1
    END SUBROUTINE VertexToElementList

    ! k-way merge with an array
    SUBROUTINE kWayMergeArray(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, map)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      INTEGER :: map(:)

      INTEGER :: i, j, k, vindi

      ! Merge nv lists using a map (i.e. an array)
      nn = 1
      DO i=1,nv
        DO j=ptrli(i), ptrti(i)-1
          vindi = vind(j)
          ! Put element to map if it is not already there
          IF (map(vindi)==0 .AND. vindi /= node) THEN
            neighind(nn)=vindi
            ! Increase counter
            map(vindi)=1
            nn=nn+1
          END IF
        END DO
      END DO
      nn=nn-1

      ! Clear map
      DO i=1,nn
        map(neighind(i)) = 0
      END DO
    END SUBROUTINE kWayMergeArray

    ! k-way merge with an actual heap
    SUBROUTINE kWayMergeHeap(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, heap)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      TYPE(IntTuple_t) :: heap(:)

      TYPE(IntTuple_t) :: tmp
      INTEGER :: ii, l, r, mind, ll, tmpval, tmpind

      ! Local variables
      INTEGER :: i, e, nzheap, vindi, lindi, pind

      ! Put elements to heap
      nzheap = 0
      DO i=1,nv
        IF (ptrli(i)<ptrti(i)) THEN
          heap(i) % i1 = vind(ptrli(i))
          heap(i) % i2= i
          ptrli(i) = ptrli(i)+1
          nzheap = nzheap+1
        END IF
      END DO

      ! Build heap
      DO ii=(nzheap/2), 1, -1
        i = ii
        ! CALL BinaryHeapHeapify(heap, nzheap, i)
        DO 
          ! Find index of the minimum element
          IF (2*i<=nzheap) THEN
            IF (heap(2*i) % i1 < heap(i) % i1) THEN
              mind = 2*i
            ELSE
              mind = i
            END IF
            IF (2*i+1<=nzheap) THEN
              IF (heap(2*i+1) % i1 < heap(mind) % i1) mind = 2*i+1
            END IF
          ELSE
            mind = i
          END IF

          IF (mind == i) EXIT

          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO
      END DO

      pind = -1
      nn = 1
      DO e=1,te
        ! Pick the first element from heap
        vindi = heap(1) % i1
        lindi = heap(1) % i2

        ! Remove duplicates
        IF (vindi /= pind .AND. vindi /= node) THEN
          neighind(nn) = vindi
          pind = vindi
          nn = nn+1
        END IF

        ! Add new element from list (if any)
        IF (ptrli(lindi) < ptrti(lindi)) THEN
          heap(1) % i1 = vind(ptrli(lindi))
          heap(1) % i2 = lindi
          ptrli(lindi) = ptrli(lindi)+1
        ELSE
          heap(1) % i1 = heap(nzheap) % i1
          heap(1) % i2 = heap(nzheap) % i2
          nzheap=nzheap-1
        END IF
        ! CALL BinaryHeapHeapify(heap, nzheap, 1)
        i = 1

        DO 
          ! Find the index of the minimum element
          ii = 2*i
          mind = i
          IF (ii+1<=nzheap) THEN
            ! Elements 2*i and 2*i+1 can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
            IF (heap(ii+1) % i1 < heap(mind) % i1) mind = ii+1
          ELSE IF (ii<=nzheap) THEN
            ! Element ii can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
          END IF

          IF (mind == i) EXIT

          ! Bubble down the element
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO

      END DO
      nn=nn-1
    END SUBROUTINE kWayMergeHeap

    SUBROUTINE BinaryHeapHeapify(heap, nelem, sind)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      INTEGER, INTENT(IN) :: sind

      INTEGER :: i, l, r, mind
      TYPE(IntTuple_t) :: tmp

      i = sind
      DO
        l = 2*i
        r = 2*i+1
        ! Find index of the minimum element
        mind = i
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) mind = l
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(mind) % i1) mind = r
        END IF

        IF (mind /= i) THEN
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        ELSE
          EXIT
        END IF
      END DO
    END SUBROUTINE BinaryHeapHeapify

    FUNCTION BinaryHeapIsHeap(heap, nelem) RESULT(heaporder)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      LOGICAL :: heaporder

      INTEGER :: i, l, r

      heaporder = .TRUE.

      DO i=(nelem/2), 1, -1
        l = 2*i
        r = 2*i+1
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'left: ', l, i
            EXIT
          END IF
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'right: ', r, i
            EXIT
          END IF
        END IF
      END DO
    END FUNCTION BinaryHeapIsHeap

  END SUBROUTINE ElmerMeshToDualGraph

  SUBROUTINE Graph_Deallocate(Graph)
    IMPLICIT NONE
    TYPE(Graph_t) :: Graph

    DEALLOCATE(Graph % ptr)
    DEALLOCATE(Graph % ind)
    Graph % n = 0
  END SUBROUTINE Graph_Deallocate

  SUBROUTINE ElmerGraphColour(Graph, Colouring, ConsistentColours)
    IMPLICIT NONE

    TYPE(Graph_t), INTENT(IN) :: Graph
    TYPE(Graphcolour_t) :: Colouring
    LOGICAL, OPTIONAL :: ConsistentColours

    INTEGER, ALLOCATABLE :: uncolored(:)
    INTEGER, ALLOCATABLE :: fc(:), ucptr(:), rc(:), rcnew(:)

    INTEGER :: nc, dualmaxdeg, i, v, w, uci, wci, vli, vti, vcol, wcol, &
            nrc, nunc, nthr, TID, allocstat, gn
    INTEGER, ALLOCATABLE :: colours(:)
    INTEGER, PARAMETER :: VERTEX_PER_THREAD = 100
    LOGICAL :: consistent

    ! Iterative parallel greedy algorithm (Alg 2.) from 
    ! U. V. Catalyurek, J. Feo, A.H. Gebremedhin, M. Halappanavar, A. Pothen. 
    ! "Graph coloring algorithms for multi-core and massively multithreaded systems".
    ! Parallel computing, 38, 2012, pp. 576--594. 

    ! Initialize number of colours, maximum degree of graph and number of 
    ! uncolored vertices
    nc = 0
    dualmaxdeg = 0
    gn = Graph % n
    nunc = gn

    ! Check if a reproducible colouring is being requested
    consistent = .FALSE.
    IF (PRESENT(ConsistentColours)) consistent = ConsistentColours

    ! Get maximum vertex degree of the given graph
    !$OMP PARALLEL DO SHARED(Graph) &
    !$OMP PRIVATE(v) REDUCTION(max:dualmaxdeg) DEFAULT(NONE)
    DO v=1,Graph % n
      dualmaxdeg = MAX(dualmaxdeg, Graph % ptr(v+1)- Graph % ptr(v))
    END DO
    !$OMP END PARALLEL DO

    nthr = 1
    ! Ensure that each vertex has at most one thread attached to it
    !$ IF (.NOT. consistent) nthr = MIN(omp_get_max_threads(), gn)

    ! Allocate memory for colours of vertices and thread colour pointers
    ALLOCATE(colours(gn), uncolored(gn), ucptr(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate colour maps!')

    !$OMP PARALLEL SHARED(gn, dualmaxdeg, Graph, colours, nunc, &
    !$OMP                 uncolored, ucptr, nthr) &
    !$OMP PRIVATE(uci, vli, vti, v, w, wci, vcol, wcol, fc, nrc, rc, rcnew, &
    !$OMP         allocstat, TID) &
    !$OMP REDUCTION(max:nc) DEFAULT(NONE) NUM_THREADS(nthr)

    TID=1
    !$ TID=OMP_GET_THREAD_NUM()+1

    ! Greedy algorithm colours a given graph with at 
    ! most max_{v\in V} deg(v)+1 colours
    ALLOCATE(fc(dualmaxdeg+1), rc((gn/nthr)+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate local workspace!')
    ! Initialize forbidden colour array (local to thread)
    fc = 0

    ! Initialize colours and uncolored entries
    !$OMP DO 
    DO v=1,gn
      colours(v)=0
      ! U <- V
      uncolored(v)=v
    END DO
    !$OMP END DO

    DO
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1

        ! For each w\in adj(v) do
        DO w=vli, vti
          ! fc[colour[w]]<-v
          !$OMP ATOMIC READ
          wcol = colours(Graph % ind(w))
          IF (wcol /= 0) fc(wcol) = v
        END DO

        ! Find smallest permissible colour for vertex
        ! c <- min\{i>0: fc[i]/=v \}
        DO i=1,dualmaxdeg+1
          IF (fc(i) /= v) THEN
            !$OMP ATOMIC WRITE 
            colours(v) = i
            ! Maintain maximum colour
            nc = MAX(nc, i)
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO

      nrc = 0
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1
        vcol = colours(v)

        ! Make sure that recolour array has enough storage for 
        ! the worst case (all elements need to be added)
        IF (SIZE(rc)<nrc+(vti-vli)+1) THEN
          ALLOCATE(rcnew(MAX(SIZE(rc)*2, nrc+(vti-vli)+1)), STAT=allocstat)
          IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
                  'Unable to allocate local workspace!')
          rcnew(1:nrc)=rc(1:nrc)
          DEALLOCATE(rc)
          CALL MOVE_ALLOC(rcnew, rc)
        END IF

        ! For each w\in adj(v) do
        DO wci=vli,vti
          w = Graph % ind(wci)
          IF (colours(w)==vcol .AND. v>w) THEN
            ! R <- R\bigcup {v} (thread local)
            nrc = nrc + 1
            rc(nrc)=v
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO NOWAIT

      ucptr(TID)=nrc
      !$OMP BARRIER

      !$OMP SINGLE
      CALL ComputeCRSIndexes(nthr, ucptr)
      nunc = ucptr(nthr+1)-1
      !$OMP END SINGLE

      ! U <- R
      uncolored(ucptr(TID):ucptr(TID+1)-1)=rc(1:nrc)
      !$OMP BARRIER

      ! Colour the remaining vertices sequentially if the 
      ! size of the set of uncoloured vertices is small enough
      IF (nunc < nthr*VERTEX_PER_THREAD) THEN
        !$OMP SINGLE
        DO uci=1,nunc
          v = uncolored(uci)
          vli = Graph % ptr(v)
          vti = Graph % ptr(v+1)-1

          ! For each w\in adj(v) do
          DO w=vli, vti
            ! fc[colour[w]]<-v
            wcol = colours(Graph % ind(w))
            IF (wcol /= 0) fc(wcol) = v
          END DO

          ! Find smallest permissible colour for vertex
          ! c <- min\{i>0: fc[i]/=v \}
          DO i=1,dualmaxdeg+1
            IF (fc(i) /= v) THEN
              ! Single thread, no collisions possible 
              colours(v) = i
              ! Maintain maximum colour
              nc = MAX(nc, i)
              EXIT
            END IF
          END DO
        END DO
        !$OMP END SINGLE NOWAIT

        EXIT
      END IF

    END DO

    ! Deallocate thread local storage
    DEALLOCATE(fc, rc)
    !$OMP END PARALLEL

    DEALLOCATE(uncolored, ucptr)

    ! Set up colouring data structure
    Colouring % nc = nc
    CALL MOVE_ALLOC(colours, Colouring % colours)
  END SUBROUTINE ElmerGraphColour

  SUBROUTINE Colouring_Deallocate(Colours)
    IMPLICIT NONE
    TYPE(GraphColour_t) :: Colours

    DEALLOCATE(Colours % colours)
    Colours % nc = 0
  END SUBROUTINE Colouring_Deallocate

  SUBROUTINE ElmerColouringToGraph(Colours, PackedList)
    IMPLICIT NONE

    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(Graph_t) :: PackedList

    INTEGER, ALLOCATABLE :: cptr(:), cind(:)

    INTEGER :: nc, c, i, n, allocstat

    nc = Colours % nc
    n = size(Colours % colours)
    ALLOCATE(cptr(nc+1), cind(n), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerGatherColourLists','Memory allocation failed.')
    cptr = 0
    ! Count number of elements in each colour
    DO i=1,n
      cptr(Colours % colours(i))=cptr(Colours % colours(i))+1
    END DO

    CALL ComputeCRSIndexes(nc, cptr)

    DO i=1,n
      c=Colours % colours(i)
      cind(cptr(c))=i
      cptr(c)=cptr(c)+1
    END DO

    DO i=nc,2,-1
      cptr(i)=cptr(i-1)
    END DO
    cptr(1)=1

    ! Set up graph data structure
    PackedList % n = nc
    CALL MOVE_ALLOC(cptr, PackedList % ptr)
    CALL MOVE_ALLOC(cind, PackedList % ind)
  END SUBROUTINE ElmerColouringToGraph

  ! Routine constructs colouring for boundary mesh based on colours of main mesh
  SUBROUTINE ElmerBoundaryGraphColour(Mesh, Colours, BoundaryColours)
    IMPLICIT NONE

    TYPE(Mesh_t), INTENT(IN) :: Mesh
    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(GraphColour_t) :: BoundaryColours

    TYPE(Element_t), POINTER :: Element
    INTEGER :: elem, nelem, nbelem, astat, lcolour, rcolour, nbc
    INTEGER, ALLOCATABLE :: bcolours(:)

    nelem = Mesh % NumberOfBulkElements
    nbelem = Mesh % NumberOfBoundaryElements

    ! Allocate boundary colouring
    ALLOCATE(bcolours(nbelem), STAT=astat)
    IF (astat /= 0) THEN
       CALL Fatal('ElmerBoundaryGraphColour','Unable to allocate boundary colouring')
    END IF
    
    nbc = 0
    ! Loop over boundary mesh
    !$OMP PARALLEL DO &
    !$OMP SHARED(Mesh, nelem, nbelem, Colours, bcolours) &
    !$OMP PRIVATE(Element, lcolour, rcolour) &
    !$OMP REDUCTION(max:nbc) &
    !$OMP DEFAULT(NONE)
    DO elem=1,nbelem       
       Element => Mesh % Elements(nelem+elem)

       ! Try to find colour for boundary element based on left / right parent
       lcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          lcolour = Colours % colours(Element % BoundaryInfo % Left % ElementIndex)
       END IF
       rcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          rcolour = Colours % colours(Element % BoundaryInfo % Right % ElementIndex)
       END IF

       ! Sanity check for debug
       IF (ASSOCIATED(Element % BoundaryInfo % Left) .AND. & 
          ASSOCIATED(Element % BoundaryInfo % Right) .AND. &
            lcolour /= rcolour) THEN
         CALL Warn('ElmerBoundaryGraphColour','Inconsistent colours for boundary element: ' &
               // TRIM(i2s(elem)) // "=>" &
               // TRIM(i2s(lcolour))// " | "//TRIM(i2s(rcolour)))
         WRITE (*,*) Element % BoundaryInfo % Left % ElementIndex, Element % BoundaryInfo % Right % ElementIndex
       END IF

       bcolours(elem)=MAX(lcolour,rcolour)
       nbc=MAX(nbc,bcolours(elem))
    END DO
    !$OMP END PARALLEL DO

    ! Set up colouring data structure
    BoundaryColours % nc = nbc
    CALL MOVE_ALLOC(bcolours, BoundaryColours % colours)
  END SUBROUTINE ElmerBoundaryGraphColour
  
  ! Given CRS indices, referenced indirectly from graph, 
  ! evenly load balance the work among the nthr threads
  SUBROUTINE ThreadLoadBalanceElementNeighbour(nthr, gn, gptr, gind, &
          rptr, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER :: gptr(:), gind(:), rptr(:)
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, j, k, wrk, gwrk, thrwrk, allocstat

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadLoadBalanceElementNeighbour', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Compute total global work
    gwrk = 0
    DO i=1,gn
      DO j=gptr(i),gptr(i+1)-1
        gwrk = gwrk + (rptr(gind(j)+1)-rptr(gind(j)))
      END DO
    END DO

    ! Amount of work per thread
    thrwrk = CEILING(REAL(gwrk,dp) / nthr)

    ! Find rows for each thread to compute
    blkleads(1)=1
    DO i=1,nthr
      wrk = 0
      ! Acquire enough work for thread i
      DO j=blkleads(i),gn
        DO k=gptr(j),gptr(j+1)-1
          wrk = wrk + (rptr(gind(j)+1)-rptr(gind(j)))
        END DO
        IF (wrk >= thrwrk) EXIT
      END DO

      blkleads(i+1)=j+1
      ! Check if we have run out of rows
      IF (j+1>gn) EXIT
    END DO
    ! Reset number of rows (may be less than or equal to original number)
    nthr = i
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadLoadBalanceElementNeighbour

  SUBROUTINE ThreadStaticWorkShare(nthr, gn, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, rem, thrwrk, allocstat
    INTEGER :: totelem

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadStaticWorkShare', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Assuming even distribution of nodes / element, 
    ! distribute rows for each thread to compute 
    blkleads(1)=1
    thrwrk = gn / nthr
    rem = gn-nthr*thrwrk
    ! totelem = 0
    DO i=1,nthr-1
      IF (i<rem) THEN
        blkleads(i+1)=blkleads(i)+thrwrk+1
      ELSE
        blkleads(i+1)=blkleads(i)+thrwrk
      END IF
    END DO
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadStaticWorkShare

  ! Given row counts, in-place compute CRS indices to data
  SUBROUTINE ComputeCRSIndexes(n, arr)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER :: arr(:)

    INTEGER :: i, indi, indip

    indi = arr(1)
    arr(1)=1
    DO i=1,n-1
      indip=arr(i+1)
      arr(i+1)=arr(i)+indi
      indi=indip
    END DO
    arr(n+1)=arr(n)+indi
  END SUBROUTINE ComputeCRSIndexes

  !> Calculate body average for a discontinuous galerkin field.
  !> The intended use is in conjunction of saving the results. 
  !> This tampers the field and therefore may have unwanted side effects
  !> if the solution is to be used for something else too.
  !-------------------------------------------------------------------
  SUBROUTINE CalculateBodyAverage( Mesh, Var, BodySum )

    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: BodySum

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: BodyAverage(:)
    INTEGER, ALLOCATABLE :: BodyCount(:)
    INTEGER :: n,i,j,k,l,nodeind,dgind, Nneighbours
    REAL(KIND=dp) :: AveHits
    LOGICAL, ALLOCATABLE :: IsNeighbour(:)
    LOGICAL :: Parallel

    
    IF(.NOT. ASSOCIATED(var)) RETURN
    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) RETURN

    IF( Var % DgAveraged ) THEN
      CALL Info('CalculateBodyAverage','Nodal average already computed for: '&
          //TRIM(Var % Name), Level=15)
      RETURN
    END IF
    
    IF( BodySum ) THEN
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal sum for: '&
          //TRIM(Var % Name), Level=8)
    ELSE
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal average for: '&
          //TRIM(Var % Name), Level=8)
    END IF

    Parallel = (ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 
    
    
    n = Mesh % NumberOfNodes
    ALLOCATE( BodyCount(n), BodyAverage(n), IsNeighbour(Parenv % PEs) )
  
    
    DO i=1,CurrentModel % NumberOfBodies

      DO k=1,Var % Dofs
        BodyCount = 0
        BodyAverage = 0.0_dp

        DO j=1,Mesh % NumberOfBulkElements 
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE
          DO l = 1, Element % TYPE % NumberOfNodes
            nodeind = Element % NodeIndexes(l)
            dgind = Var % Perm(Element % DGIndexes(l) )
            IF( dgind > 0 ) THEN
              BodyAverage( nodeind ) = BodyAverage( nodeind ) + &
                  Var % Values( Var % DOFs*( dgind-1)+k )
              BodyCount( nodeind ) = BodyCount( nodeind ) + 1 
            END IF
          END DO
        END DO

        IF( k == 1 ) THEN
          ! This is just low priority info on the averaging
          IF( InfoActive(20) ) THEN
            j = COUNT(BodyCount > 0) 
            IF( j > 0 ) THEN
              AveHits = 1.0_dp * SUM( BodyCount ) / j
            ELSE
              AveHits = 0.0_dp
            END IF
            WRITE(Message,'(A,ES12.3)') 'In body '//TRIM(I2S(i))//' average hit count is: ',AveHits
            CALL Info('CalculateBodyAverage',Message) 
            WRITE(Message,'(A,2I0)') 'In body '//TRIM(I2S(i))//' hit count range is: ',&
                MINVAL(BodyCount,BodyCount>0), MAXVAL(BodyCount)
            CALL Info('CalculateBodyAverage',Message) 
          END IF
        END IF
          
        IF( Parallel ) THEN
          Nneighbours = MeshNeighbours(Mesh, IsNeighbour)
          CALL SendInterface(); CALL RecvInterface()
        END IF

        ! Do not average weighted quantities (like nodal forces) - they should only be summed.
        ! But do average all other quantities. 
        IF( .NOT. BodySum ) THEN
          DO j=1,n
            IF( BodyCount(j) > 0 ) BodyAverage(j) = BodyAverage(j) / BodyCount(j)
          END DO
        END IF

        ! Now copy the average values to the DG field
        DO j=1,Mesh % NumberOfBulkElements 
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE
          DO l = 1, Element % TYPE % NumberOfNodes
            nodeind = Element % NodeIndexes(l)
            dgind = Var % Perm(Element % DGIndexes(l) )
            IF( dgind > 0 ) THEN
              Var % Values( Var % DOFs*( dgind-1)+k ) = BodyAverage( nodeind ) 
            END IF
          END DO
        END DO
      END DO
    END DO

    Var % DgAveraged = .TRUE.
    
CONTAINS

     SUBROUTINE SendInterface()
       TYPE buf_t
         REAL(KIND=dp), ALLOCATABLE :: dval(:)
         INTEGER, ALLOCATABLE :: gdof(:), ival(:)
       END TYPE buf_t

       INTEGER, ALLOCATABLE :: cnt(:)
       TYPE(buf_t), ALLOCATABLE :: buf(:)

       INTEGER :: i,j,k,ierr

       ALLOCATE(cnt(ParEnv % PEs), buf(ParEnv % PEs))

       cnt = 0
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT.Mesh % ParallelInfo % NodeInterface(i)) CYCLE
         IF(BodyCount(i) <= 0 ) CYCLE

         DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
           k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)+1
           cnt(k) = cnt(k) + 1
         END DO
       END DO

       DO i=1,ParEnv % PEs
         ALLOCATE(buf(i) % gdof(cnt(i)), buf(i) % ival(cnt(i)), buf(i) % dval(cnt(i)))
       END DO

       cnt = 0
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT.Mesh % ParallelInfo % NodeInterface(i)) CYCLE
         IF(BodyCount(i) <= 0 ) CYCLE

         DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
           k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)+1
           cnt(k) = cnt(k) + 1
           buf(k) % gdof(cnt(k)) = Mesh % ParallelInfo % GlobalDOFs(i)
           buf(k) % ival(cnt(k)) = BodyCount(i)
           buf(k) % dval(cnt(k)) = BodyAverage(i)
         END DO
       END DO

       DO i=1,ParEnv % PEs
         IF(.NOT. isNeighbour(i)) CYCLE

         CALL MPI_BSEND( cnt(i),1,MPI_INTEGER,i-1,1310,ELMER_COMM_WORLD,ierr )
         IF(cnt(i)>0) THEN
           CALL MPI_BSEND( buf(i) % gdof,cnt(i),MPI_INTEGER,i-1,1311,ELMER_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % ival,cnt(i),MPI_INTEGER,i-1,1312,ELMER_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % dval,cnt(i),MPI_DOUBLE_PRECISION,i-1,1313,ELMER_COMM_WORLD,ierr )
         END IF
       END DO
     END SUBROUTINE SendInterface


     SUBROUTINE RecvInterface()
       INTEGER, ALLOCATABLE :: gdof(:), ival(:)
       REAL(KIND=dp), ALLOCATABLE :: dval(:)
       INTEGER :: i,j,k,ierr, cnt, status(MPI_STATUS_SIZE)

       DO i=1,ParEnv % PEs

         IF(.NOT.isNeighbour(i)) CYCLE

         CALL MPI_RECV( cnt,1,MPI_INTEGER,i-1,1310,ELMER_COMM_WORLD,status,ierr )
         IF(cnt>0) THEN
           ALLOCATE( gdof(cnt), ival(cnt), dval(cnt) )
           CALL MPI_RECV( gdof,cnt,MPI_INTEGER,i-1,1311,ELMER_COMM_WORLD,status,ierr )
           CALL MPI_RECV( ival,cnt,MPI_INTEGER,i-1,1312,ELMER_COMM_WORLD,status,ierr )
           CALL MPI_RECV( dval,cnt,MPI_DOUBLE_PRECISION,i-1,1313,ELMER_COMM_WORLD,status,ierr )

           DO j=1,cnt
             k = SearchNode(Mesh % ParallelInfo, gdof(j))
             IF (k>0) THEN
               BodyCount(k) = BodyCount(k) + ival(j)
               BodyAverage(k) = BodyAverage(k)  + dval(j)
             END IF
           END DO 
           DEALLOCATE( gdof, ival, dval )
         END IF
       END DO
       CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
     END SUBROUTINE RecvInterface

  END SUBROUTINE CalculateBodyAverage



  !> Given an elemental DG field create a minimal reduced set of it that maintains
  !> the necessary continuities. The continuities may be requested between bodies
  !> or materials. Optionally the user may give a boundary mask which defines the 
  !> potential discontinuous nodes that may be greedy or not. 
  !-------------------------------------------------------------------------------
  FUNCTION MinimalElementalSet( Mesh, JumpMode, VarPerm, BcFlag, &
      NonGreedy ) RESULT ( SetPerm )

    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: JumpMode
    INTEGER, POINTER, OPTIONAL :: VarPerm(:)
    CHARACTER(LEN=*), OPTIONAL :: BcFlag
    LOGICAL, OPTIONAL :: NonGreedy
    INTEGER, POINTER :: SetPerm(:)

    TYPE(Element_t), POINTER :: Element, Left, Right
    INTEGER :: n,i,j,k,l,bc_id,mat_id,body_id,NoElimNodes,nodeind,JumpModeIndx,&
        LeftI,RightI,NumberOfBlocks
    LOGICAL, ALLOCATABLE :: JumpNodes(:)
    INTEGER, ALLOCATABLE :: NodeVisited(:)
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Found
    

    CALL Info('MinimalElementalSet','Creating discontinuous subset from DG field',Level=5)

    ! Calculate size of permutation vector
    ALLOCATE( NodeVisited( Mesh % NumberOfNodes ) )
    NodeVisited = 0

    NULLIFY( SetPerm ) 
    k = 0
    DO i=1,Mesh % NumberOfBulkElements         
      Element => Mesh % Elements(i)
      k = k + Element % TYPE % NumberOfNodes
    END DO
    CALL Info('MinimalElementalSet','Maximum number of dofs in DG: '//TRIM(I2S(k)),Level=12)
    ALLOCATE( SetPerm(k) )
    SetPerm = 0
    l = 0
    NoElimNodes = 0

    CALL Info('MinimalElementalSet','Reducing elemental discontinuity with mode: '//TRIM(JumpMode),Level=7)

    SELECT CASE ( JumpMode )

    CASE('db') ! discontinuous bodies
      NumberOfBlocks = CurrentModel % NumberOfBodies
      JumpModeIndx = 1

    CASE('dm') ! discontinuous materials
      NumberOfBlocks = CurrentModel % NumberOfMaterials
      JumpModeIndx = 2

    CASE DEFAULT
      CALL Fatal('MinimalElementalSet','Unknown JumpMode: '//TRIM(JumpMode))

    END SELECT
  

    IF( PRESENT( BcFlag ) ) THEN
      ALLOCATE( JumpNodes( Mesh % NumberOfNodes ) )
    END IF

    
    DO i=1,NumberOfBlocks
      
      ! Before the 1st block no numbers have been given.
      ! Also if we want discontinuous blocks on all sides initialize the whole list to zero. 
      IF( i == 1 .OR. .NOT. PRESENT( BcFlag ) ) THEN
        NodeVisited = 0

      ELSE
        ! Vector indicating the disontinuous nodes
        ! If this is not given all interface nodes are potentially discontinuous
        JumpNodes = .FALSE.
        
        DO j=Mesh % NumberOfBulkElements + 1, &
            Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(j)

          DO bc_id=1,CurrentModel % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
          END DO
          IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE
          IF( .NOT. ListCheckPresent( CurrentModel % BCs(bc_id) % Values, BcFlag ) ) CYCLE

          Left => Element % BoundaryInfo % Left
          Right => Element % BoundaryInfo % Right
          IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) CYCLE

          IF( JumpModeIndx == 1 ) THEN
            LeftI = Left % BodyId
            RightI = Right % BodyId
          ELSE
            LeftI = ListGetInteger( CurrentModel % Bodies(Left % BodyId) % Values,'Material',Found)
            RightI = ListGetInteger( CurrentModel % Bodies(Right % BodyId) % Values,'Material',Found)
          END IF

          IF( LeftI /= i .AND. RightI /= i ) CYCLE
          JumpNodes( Element % NodeIndexes ) = .TRUE.
        END DO

        IF( PRESENT( NonGreedy ) ) THEN
          IF( NonGreedy ) THEN        
            DO j=Mesh % NumberOfBulkElements + 1, &
                Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
              Element => Mesh % Elements(j)

              DO bc_id=1,CurrentModel % NumberOfBCs
                IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
              END DO
              IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE

              IF( ListCheckPresent( CurrentModel % BCs(bc_id) % Values, BcFlag ) ) CYCLE

              Left => Element % BoundaryInfo % Left
              Right => Element % BoundaryInfo % Right

              ! External BCs don't have a concept of jump, so no need to treat them
              IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) CYCLE

              JumpNodes( Element % NodeIndexes ) = .FALSE.
            END DO
          END IF
        END IF

        ! Initialize new potential nodes for the block where we found discontinuity
        WHERE( JumpNodes ) NodeVisited = 0
      END IF


      ! Now do the real thing. 
      ! Add new dofs such that minimal discontinuity is maintained 
      DO j=1,Mesh % NumberOfBulkElements         
        Element => Mesh % Elements(j)

        Body_Id = Element % BodyId 
        IF( JumpModeIndx == 1 ) THEN
          IF( Body_id /= i ) CYCLE
        ELSE
          Mat_Id = ListGetInteger( CurrentModel % Bodies(Body_Id) % Values,'Material',Found)
          IF( Mat_Id /= i ) CYCLE
        END IF

        NodeIndexes => Element % NodeIndexes
        
        DO k=1,Element % TYPE % NumberOfNodes         
          nodeind = NodeIndexes(k)
          IF( PRESENT( VarPerm ) ) THEN
            IF( VarPerm( nodeind ) == 0 ) CYCLE
          END IF
          IF( NodeVisited( nodeind ) > 0 ) THEN
            SetPerm( Element % DGIndexes(k) ) = NodeVisited( nodeind )
            NoElimNodes = NoElimNodes + 1
          ELSE
            l = l + 1
            NodeVisited(nodeind) = l
            SetPerm( Element % DGIndexes(k) ) = l
          END IF
        END DO
      END DO
    END DO

    CALL Info('MinimalElementalSet','Independent dofs in elemental field: '//TRIM(I2S(l)),Level=7)
    CALL Info('MinimalElementalSet','Redundant dofs in elemental field: '//TRIM(I2S(NoElimNodes)),Level=7)     

  END FUNCTION MinimalElementalSet


  !> Calculate the reduced DG field given the reduction permutation.
  !> The permutation must be predefined. This may be called repeatedly
  !> for different variables. Optionally one may take average, or 
  !> a plain sum over the shared nodes. 
  !-------------------------------------------------------------------
  SUBROUTINE ReduceElementalVar( Mesh, Var, SetPerm, TakeAverage )

    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: SetPerm(:)
    LOGICAL :: TakeAverage

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: SetSum(:)
    INTEGER, ALLOCATABLE :: SetCount(:)
    INTEGER :: dof,n,m,i,j,k,l,nodeind,dgind
    REAL(KIND=dp) :: AveHits

    IF(.NOT. ASSOCIATED(var)) THEN
      CALL Warn('ReduceElementalVar','Variable not associated!')
      RETURN
    END IF

    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) THEN
      CALL Warn('ReduceElementalVar','Var % Perm too small!')
      RETURN
    END IF

    IF( TakeAverage ) THEN
      CALL Info('ReduceElementalVar','Calculating reduced set average for: '&
          //TRIM(Var % Name), Level=7)
    ELSE
      CALL Info('ReduceElementalVar','Calculating reduced set sum for: '&
          //TRIM(Var % Name), Level=7)
    END IF

    n = Mesh % NumberOfNodes

    m = MAXVAL( SetPerm )
    ALLOCATE( SetCount(m), SetSum(m) )
    SetCount = 0
    SetSum = 0.0_dp

    ! Take the sum to nodes, and calculate average if requested
    DO dof=1,Var % Dofs
      SetCount = 0
      SetSum = 0.0_dp

      DO i=1,SIZE(SetPerm)
        j = SetPerm(i)
        l = Var % Perm(i)
        SetSum(j) = SetSum(j) + Var % Values( Var % DOFs * (l-1) + dof )
        SetCount(j) = SetCount(j) + 1
      END DO
        
      m = SUM( SetCount ) 
      IF( m == 0 ) RETURN

      IF( TakeAverage ) THEN
        WHERE( SetCount > 0 ) SetSum = SetSum / SetCount
      END IF

      IF( dof == 1 ) THEN
        AveHits = 1.0_dp * SUM( SetCount ) / COUNT( SetCount > 0 )
        WRITE(Message,'(A,ES15.4)') 'Average number of hits: ',AveHits
        CALL Info('ReduceElementalVar',Message,Level=10)
      END IF

      ! Copy the reduced set back to the original elemental field
      DO i=1,SIZE(SetPerm)
        j = SetPerm(i)
        l = Var % Perm(i)
        Var % Values( Var % DOFs * (l-1) + dof ) = SetSum(j)
      END DO
    END DO

  END SUBROUTINE ReduceElementalVar


  !> Given a elemental DG field and a reduction permutation compute the 
  !> body specific lumped sum. The DG field may be either original one
  !> or already summed up. In the latter case only one incident of the 
  !> redundant nodes is set.
  !---------------------------------------------------------------------
  SUBROUTINE LumpedElementalVar( Mesh, Var, SetPerm, AlreadySummed )
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: SetPerm(:)
    LOGICAL :: AlreadySummed

    TYPE(Element_t), POINTER :: Element
    LOGICAL, ALLOCATABLE :: NodeVisited(:)
    INTEGER :: dof,n,m,i,j,k,l,nodeind,dgind
    REAL(KIND=dp), ALLOCATABLE :: BodySum(:)

    IF(.NOT. ASSOCIATED(var)) RETURN
    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) RETURN

    CALL Info('LumpedElementalVar','Calculating lumped sum for: '&
        //TRIM(Var % Name), Level=8)

    n = Mesh % NumberOfNodes

    m = MAXVAL( SetPerm )
    IF( AlreadySummed ) THEN
      ALLOCATE( NodeVisited(m) )
    END IF
    ALLOCATE( BodySum( CurrentModel % NumberOfBodies ) )

    ! Take the sum to nodes, and calculate average if requested
    DO dof=1,Var % Dofs

      BodySum = 0.0_dp

      DO i=1,CurrentModel % NumberOfBodies

        IF( AlreadySummed ) THEN
          NodeVisited = .FALSE.
        END IF

        DO j=1,Mesh % NumberOfBulkElements         
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE

          DO k=1,Element % TYPE % NumberOfNodes         
            dgind = Element % DGIndexes(k)
            l = SetPerm(dgind)
            IF( l == 0 ) CYCLE

            IF( AlreadySummed ) THEN
              IF( NodeVisited(l) ) CYCLE           
              NodeVisited(l) = .TRUE.
            END IF

            BodySum(i) = BodySum(i) + &
                Var % Values( Var % Dofs * ( Var % Perm( dgind )-1) + dof )
          END DO
        END DO
      END DO

      IF( Var % Dofs > 1 ) THEN
        CALL Info('LumpedElementalVar','Lumped sum for component: '//TRIM(I2S(dof)),Level=6)
      END IF
      DO i=1,CurrentModel % NumberOfBodies
        WRITE(Message,'(A,ES15.4)') 'Body '//TRIM(I2S(i))//' sum:',BodySum(i)
        CALL Info('LumpedElementalVar',Message,Level=10)
      END DO

    END DO

    DEALLOCATE( NodeVisited, BodySum )

  END SUBROUTINE LumpedElementalVar



!------------------------------------------------------------------------------
  SUBROUTINE SaveParallelInfo( Solver )
!------------------------------------------------------------------------------
   TYPE( Solver_t ), POINTER  :: Solver
!------------------------------------------------------------------------------    
   TYPE(ParallelInfo_t), POINTER :: ParInfo=>NULL()
   TYPE(ValueList_t), POINTER :: Params
   CHARACTER(LEN=MAX_NAME_LEN) :: dumpfile
   INTEGER :: i,j,k,n,maxnei
   LOGICAL :: Found, MeshMode, MatrixMode
   CHARACTER(*), PARAMETER :: Caller = "SaveParallelInfo"
   TYPE(Nodes_t), POINTER :: Nodes
   
   Params => Solver % Values 

   MeshMode = ListGetLogical( Params,'Save Parallel Matrix Info',Found ) 
   MatrixMode = ListGetLogical( Params,'Save Parallel Mesh Info',Found ) 

   IF( .NOT. ( MeshMode .OR. MatrixMode ) ) RETURN

10 IF( MeshMode ) THEN
     CALL Info(Caller,'Saving parallel mesh info',Level=8 ) 
   ELSE
     CALL Info(Caller,'Saving parallel matrix info',Level=8 ) 
   END IF

   IF( MeshMode ) THEN
     ParInfo => Solver % Mesh % ParallelInfo
     Nodes => Solver % Mesh % Nodes
     dumpfile = 'parinfo_mesh.dat'
   ELSE
     ParInfo => Solver % Matrix % ParallelInfo
     dumpfile = 'parinfo_mat.dat'      
   END IF

   IF( .NOT. ASSOCIATED( ParInfo ) ) THEN
     CALL Warn(Caller,'Parallel info not associated!')
     RETURN
   END IF

   n = SIZE( ParInfo % GlobalDOFs )
   IF( n <= 0 ) THEN
     CALL Warn(Caller,'Parallel info size is invalid!')
     RETURN
   END IF

   ! memorize the maximum number of parallel neighbours
   maxnei = 0
   IF( ASSOCIATED( ParInfo % NeighbourList ) ) THEN
     DO i=1,n
       IF( ASSOCIATED( ParInfo % NeighbourList(i) % Neighbours ) ) THEN
         j = SIZE( ParInfo % NeighbourList(i) % Neighbours )
         maxnei = MAX( j, maxnei ) 
       END IF
     END DO
   END IF
   CALL Info(Caller,'Maximum number of parallel neighbours:'//TRIM(I2S(maxnei)))

   IF(ParEnv % PEs > 1) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))      
   CALL Info(Caller,'Saving parallel info to: '//TRIM(dumpfile),Level=8)

   OPEN(1,FILE=dumpfile, STATUS='Unknown')  
   DO i=1,n
     j = ParInfo % GlobalDOFs(i)
     IF( ParInfo % NodeInterface(i) ) THEN
       k = 1
     ELSE
       k = 0
     END IF
     WRITE(1,'(3I6)',ADVANCE='NO') i,j,k
     IF( ASSOCIATED( ParInfo % NeighbourList(i) % Neighbours ) ) THEN
       k = SIZE( ParInfo % NeighbourList(i) % Neighbours )
     ELSE
       k = 0
     END IF
     DO j=1,k
       WRITE(1,'(I6)',ADVANCE='NO')  ParInfo % NeighbourList(i) % Neighbours(j)
     END DO
     DO j=k+1,maxnei
       WRITE(1,'(I6)',ADVANCE='NO')  -1 
     END DO
     IF( MeshMode ) THEN
       WRITE(1,'(3ES12.3)',ADVANCE='NO') &
           Nodes % x(i), Nodes % y(i), Nodes % z(i)
     END IF
     WRITE(1,'(A)') ' ' ! finish the line
   END DO
   CLOSE(1)

   ! Redo with matrix if both modes are requested
   IF( MeshMode .AND. MatrixMode ) THEN
     MeshMode = .FALSE.
     GOTO 10
   END IF
   
   CALL Info(Caller,'Finished saving parallel info',Level=10)

!------------------------------------------------------------------------------
 END SUBROUTINE SaveParallelInfo
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetLagrangeIndexes( Mesh, LagN, Element, Indexes )  RESULT(L)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: LagN
    TYPE(Element_t), OPTIONAL, TARGET :: Element
    INTEGER, OPTIONAL :: Indexes(:)
    INTEGER :: L
!------------------------------------------------------------------------------    
    TYPE(Solver_t),  POINTER :: Solver
    TYPE(Element_t), POINTER :: Parent, Edge, Face
    LOGICAL :: OrientationsMatch
    LOGICAL :: EdgesActive, FacesActive
    LOGICAL :: Visited = .FALSE.

    INTEGER, PARAMETER :: MAX_LAGRANGE_NODES = 729

    INTEGER :: EdgeMap(2), FaceMap(4)
    INTEGER :: VTKTetraFaceMap(4,3)
    INTEGER :: VTKBrickFaceMap(6,4), BrickFaceOrdering(6)
    INTEGER :: Perm(MAX_LAGRANGE_NODES), TmpInd(MAX_LAGRANGE_NODES)
    INTEGER :: i,j,m,n0,e1,e2,f
    INTEGER :: nelem, nface, nedge, elemdim, thiselem
    INTEGER :: nelem_max, nface_max, nedge_max
    INTEGER :: ElemFamily, nsize
    INTEGER :: ElemType
    CHARACTER(*), PARAMETER :: Caller = 'GetLagrangeIndexes'

    SAVE Visited, nelem_max, nface_max, nedge_max, nsize, EdgesActive, &
        FacesActive, VTKTetraFaceMap, VTKBrickFaceMap, BrickFaceOrdering
!------------------------------------------------------------------------------
    
    IF (.NOT. Visited) THEN
      Visited = .TRUE.

      ! VTK's convention:
      VTKTetraFaceMap(1,:) = (/ 1,2,4 /)
      VTKTetraFaceMap(2,:) = (/ 3,4,2 /)
      VTKTetraFaceMap(3,:) = (/ 1,4,3 /)
      VTKTetraFaceMap(4,:) = (/ 1,3,2 /)

      VTKBrickFaceMap(1,:) = (/ 1,4,8,5 /)
      VTKBrickFaceMap(2,:) = (/ 2,3,7,6 /)
      VTKBrickFaceMap(3,:) = (/ 1,2,6,5 /)
      VTKBrickFaceMap(4,:) = (/ 4,3,7,8 /)
      VTKBrickFaceMap(5,:) = (/ 1,2,3,4 /)
      VTKBrickFaceMap(6,:) = (/ 5,6,7,8 /)
      BrickFaceOrdering = (/ 6,4,3,5,1,2 /)

      nedge_max = 0
      nface_max = 0
      nelem_max = 0

      DO i=1,Mesh % NumberOfBulkElements
        ElemFamily = Mesh % Elements(i) % TYPE % ElementCode / 100
        CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
        nedge_max = MAX(nedge, nedge_max)
        nface_max = MAX(nface, nface_max)
        nelem_max = MAX(nelem, nelem_max) 
      END DO
      
      EdgesActive = ASSOCIATED(Mesh % Edges)
      FacesActive = ASSOCIATED(Mesh % Faces)
      
      IF (.NOT. EdgesActive .AND. nedge_max > 0) CALL Warn(Caller, 'Mesh edges needed but not associated')
      IF (.NOT. FacesActive .AND. nface_max > 0) CALL Warn(Caller, 'Mesh faces needed but not associated')

      nsize = Mesh % NumberOfNodes + nelem_max * Mesh % NumberOfBulkElements + &
          nface_max * Mesh % NumberOfFaces + nedge_max * Mesh % NumberOfEdges

      nsize = nsize + Mesh % NumberOfBoundaryElements * MAX(nedge_max, nface_max)
    END IF

    ! If we don't have a specific element, then only return the total number which is sufficiently large
    ! in order to index all DOFs in the Lagrange mesh. 
    IF (.NOT. PRESENT(Element)) THEN
      l = nsize
      RETURN
    END IF
        
    ! The count of corner nodes:
    l = Element % TYPE % ElementCode / 100 
    IF( l >= 5 .AND. l <= 7 ) l = l-1             

    IF (PRESENT(Indexes)) THEN
      Indexes = 0
      Indexes(1:l) = Element % NodeIndexes(1:l)
    END IF
    ! Offset
    n0 = Mesh % NumberOfNodes

    IF(l>4) THEN
      ElemDim = 3
    ELSE IF(l>2) THEN
      ElemDim = 2
    ELSE
      ElemDim = 1
    END IF

    
    ! Number the additional edge nodes
    IF (EdgesActive ) THEN
      ElemFamily = Element % TYPE % ElementCode / 100
      CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)

      ! If this is a boundary element, we need to number it just as it would if it were an edge
      ! of a bulk element. 
      IF( ElemDim == 1 .AND. ASSOCIATED(Element % BoundaryInfo) ) THEN
        thiselem = 1        
        nedge = nelem
      ELSE
        thiselem = 0
      END IF
      
      DO i=1,MAX(thiselem,Element % TYPE % NumberOfEdges)
        IF(thiselem==1) THEN
          ! We use sneaky definitions here to be able to use rest of the edge indexing code.
          ! We want to use the edge indexing that has been generated for the edges of
          ! the parent element.
          Parent => Element % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED( Parent ) ) THEN
            Parent => Element % BoundaryInfo % Right
          END IF
          IF (.NOT. ASSOCIATED(Parent)) RETURN
          Edge => Find_Edge(Mesh,Parent,Element)
          EdgeMap = [1,2]
        ELSE
          f = i
          SELECT CASE(ElemFamily)
          CASE(2)
            CALL Error(Caller, '2D element is supposed to have elemental DOFs')
          CASE(3)
            EdgeMap = GetTriangleEdgeMap(i)
          CASE(4)
            EdgeMap = GetQuadEdgeMap(i)
          CASE(5)
            EdgeMap = GetTetraEdgeMap(i)
            IF (i == 3) THEN
              e1 = EdgeMap(2)
              e2 = EdgeMap(1)
              EdgeMap(1) = e1
              EdgeMap(2) = e2
            END IF
          CASE(6)
            EdgeMap = GetPyramidEdgeMap(i)
          CASE(7)
            EdgeMap = GetWedgeEdgeMap(i)
          CASE(8)
            ! It seems that VTK cell types 72 and 12/29 are not interchangeable:
            IF (LagN > 2) THEN
              ! The following is needed for 72:
              SELECT CASE(i)
              CASE(11)
                f = 12
              CASE(12)
                f = 11
              CASE DEFAULT
                CONTINUE
              END SELECT
            END IF
            EdgeMap = GetBrickEdgeMap(f)
          END SELECT
          Edge => Mesh % Edges(Element % EdgeIndexes(f))     
        END IF
        
        e1 = Edge % NodeIndexes(1)
        e2 = Edge % NodeIndexes(2)

        IF (e2 < e1) THEN
          OrientationsMatch = e1 == Element % NodeIndexes(EdgeMap(2))
        ELSE
          OrientationsMatch = e1 == Element % NodeIndexes(EdgeMap(1))
        END IF        
        
        ! Ensure the edge DOFs are listed in the right order:
        IF (OrientationsMatch) THEN
          DO j=1,nedge
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = n0 + nedge_max*(Edge % ElementIndex-1)+j
          END DO
        ELSE
          DO j=nedge,1,-1
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = n0 + nedge_max*(Edge % ElementIndex-1)+j
          END DO
        END IF
      END DO

      ! Nothing to be done here. This was boundary element that was exhausted.
      IF(thiselem==1) RETURN
      
      n0 = n0 + Mesh % NumberOfEdges * nedge_max      
    END IF

    ! Then number the additional face nodes
    IF (FacesActive) THEN
      
      SELECT CASE(Element % TYPE % ElementCode / 100)
      CASE(3,4)
        ! For 2D element only save the face if it is a boundary!
        IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
          Parent => Element % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED( Parent ) ) THEN
            Parent => Element % BoundaryInfo % Right
          END IF
          IF (.NOT. ASSOCIATED(Parent)) RETURN
          Face => Find_Face(Mesh,Parent,Element)
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          
          IF (nelem < 1) RETURN

          IF (ElemFamily == 4) THEN
            Perm = LagrangeQuadFacePermutation(Element % NodeIndexes(1:4), LagN)
          ELSE
            Perm(1:3) = LagrangeTriFacePermutation(Element % NodeIndexes(1:3), LagN)
          END IF

          IF (PRESENT(Indexes)) THEN
            DO j=1,nelem
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF
          ! Permute to create the final list of indices:
          DO j=1,nelem
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END IF
        RETURN          

      CASE(5)
        DO i=1,Element % Type % NumberOfFaces
          !
          ! Elmer has created its face indices by using face maps different from
          ! VTK's convention. Set f so that we can assign the right global indices
          ! to the face i according to VTK's convention.
          !
          IF (i == 4) THEN
            f = 1
          ELSE
            f = i+1
          END IF
          
          Face => Mesh % Faces(Element % FaceIndexes(f))          
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          ! test:
          !m = 0
          !DO j=1,3
          !  DO k=1,3
          !    IF (Face % NodeIndexes(j) == Element % NodeIndexes(VTKTetraFaceMap(i,k))) THEN
          !      m = m + 1
          !      EXIT
          !    END IF
          !  END DO
          !END DO
          !IF (m /= 3) CALL Fatal(Caller, 'Face is not identified correctly')

          Perm(1:3) = LagrangeTriFacePermutation(Element % NodeIndexes(VTKTetraFaceMap(i,1:3)), LagN)
          
          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF

          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO

      CASE(6)
        ! The quad face:
        Face => Mesh % Faces(Element % FaceIndexes(1))
        CALL LagrangeDOFCount(4, LagN, nedge, nface, nelem)
        nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D
        FaceMap = GetPyramidFaceMap(1)

        IF (nface > 1) THEN
          CALL Fatal(Caller, 'For pyramids Lagrange Element Degree < 3 supported currently')
        END IF

        DO j=1,nface
          l = l + 1
          IF (PRESENT(Indexes)) Indexes(l) = n0 + nface_max*(Face % ElementIndex-1) + j
        END DO

        ! TO DO: Index triangular faces for degrees p > 3

      CASE(7)
        ! Triangular faces:
        DO f=1,2
          Face => Mesh % Faces(Element % FaceIndexes(f))
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          IF (nface < 1) CYCLE

          FaceMap = GetWedgeFaceMap(f)
          Perm(1:3) = LagrangeTriFacePermutation(Element % NodeIndexes(FaceMap(1:3)), LagN)
          
          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF

          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO

        ! Quad faces:
        DO f=3,5
          Face => Mesh % Faces(Element % FaceIndexes(f))          
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          IF (nface < 1) CYCLE

          FaceMap = GetWedgeFaceMap(f)
          Perm = LagrangeQuadFacePermutation(Element % NodeIndexes(FaceMap(1:4)), LagN)

          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF

          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO

      CASE(8)
        DO i=1,Element % Type % NumberOfFaces 
          f = BrickFaceOrdering(i)

          Face => Mesh % Faces(Element % FaceIndexes(f))          
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          IF (nface < 1) CYCLE

          Perm = LagrangeQuadFacePermutation(Element % NodeIndexes(VTKBrickFaceMap(i,1:4)), LagN)

          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF
          ! Permute to create the final list of indices:
          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO
      END SELECT
      
      n0 = n0 + Mesh % NumberOfFaces * nface_max
    END IF

    ! Then number the additional internal nodes (never shared)
    ElemFamily = Element % TYPE % ElementCode / 100
    CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)    
    DO j=1,nelem
      l = l + 1
      IF (PRESENT(Indexes)) Indexes(l) = n0 + nelem_max*(Element % ElementIndex-1) + j
    END DO

  CONTAINS
    ! 
    ! A subroutine for returning the maximal number of interior nodes associated with
    ! the element edges, faces and the volume in the Lagrange interpolation of degree p
    !
    SUBROUTINE LagrangeDOFCount(Family, p, nedge, nface, nelem)
      INTEGER, INTENT(IN) :: Family, p
      INTEGER, INTENT(OUT) :: nedge, nface, nelem
      
      INTEGER :: m

      m = p - 1
      nelem = 0
      nface = 0
      nedge = 0

      IF (Family == 1) RETURN
      
      SELECT CASE(Family)
      CASE(2)
        nelem = m
      CASE(3)
        nelem = m*(m-1)/2
        nedge = m
      CASE(4)
        nelem = m*m
        nedge = m
      CASE(5)
        nelem = m*(m-1)*(m-2)/6
        nface = m*(m-1)/2
        nedge = m
      CASE(6)
        nedge = m
        nface = m*m ! the maximum is determined by quad faces
        IF (p > 1) THEN
          IF (p==2) THEN
            nelem = 1
          ELSE
            CALL Fatal('LagrangeDOFCount', 'Cannot handle pyramids of degree > 2')
          END IF
        END IF
      CASE(7)
        nedge = m
        nface = m*m ! the maximum is determined by quad faces
        nelem = m*(m-1)/2*m
      CASE(8)
        nelem = m*m*m
        nface = m*m
        nedge = m
      CASE DEFAULT          
        CALL Fatal('LagrangeDOFCount', 'Unknown element family') 
      END SELECT
    END SUBROUTINE LagrangeDOFCount

    !
    ! A function to generate a permutation vector for indexing nodes on quad faces
    !
    FUNCTION LagrangeQuadFacePermutation(FaceNodes, p) RESULT(Perm)
      INTEGER, INTENT(IN) :: FaceNodes(4)
      INTEGER, INTENT(IN) :: p       ! the order of Lagrange interpolation
      INTEGER :: Perm(MAX_LAGRANGE_NODES)

      INTEGER, PARAMETER :: MAX_LAGRANGE_NODES = 729
      INTEGER :: AllIndices((p-1)**2)
      INTEGER :: i, j, n, i0, MinEntryInd(1)

      SELECT CASE(p)
      CASE(2)
        Perm = 0
        Perm(1) = 1

      CASE DEFAULT
        !
        ! We have 4 x 2 permutation patterns. Create a permutation
        ! vector to alter the default ordering in each case. The first face
        ! index is assigned to the node which is closest to the face corner A
        ! having the smallest global index. The next indices are created in 
        ! the direction of the face edge AB, with B the smallest possible 
        ! global index.
        !
        Perm = 0
        n = (p-1)**2
        DO i=1,n
          AllIndices(i) = i
        END DO

        MinEntryInd = MINLOC(FaceNodes(1:4))
        SELECT CASE(MinEntryInd(1))
        CASE(1)
          IF (FaceNodes(4) < FaceNodes(2)) THEN
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              Perm(i0+1:i0+p-1) = AllIndices(i:n:p-1)
            END DO
          ELSE
            Perm(1:n) = AllIndices(1:n)
          END IF

        CASE(2)
          IF (FaceNodes(3) < FaceNodes(1)) THEN
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(p-i+(j-1)*(p-1))
              END DO
            END DO
          ELSE
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(i0+p-j)
              END DO
            END DO
          END IF          

        CASE(3)
          IF (FaceNodes(4) < FaceNodes(2)) THEN
            DO i=1,n
              Perm(i) = AllIndices(n+1-i)
            END DO
          ELSE
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(n+1-i-(j-1)*(p-1))
              END DO
            END DO
          END IF

        CASE(4)
          IF (FaceNodes(1) < FaceNodes(3)) THEN
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(n-p+1+i-(j-1)*(p-1))
              END DO
            END DO
          ELSE
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(n-i*(p-1)+j)
              END DO
            END DO
          END IF
        END SELECT

      END SELECT
    END FUNCTION LagrangeQuadFacePermutation

    !
    ! A function to generate a permutation vector for indexing nodes on triangular faces
    !
    FUNCTION LagrangeTriFacePermutation(FaceNodes, p) RESULT(Perm)
      INTEGER, INTENT(IN) :: FaceNodes(3)
      INTEGER, INTENT(IN) :: p       ! the order of Lagrange interpolation
      INTEGER :: Perm(3)

      INTEGER :: MinEntryInd(1)

      SELECT CASE(p)
      CASE(3)
        Perm = 0
        Perm(1) = 1

      CASE(4)
        !
        ! We have 3 x 2 permutation patterns. Create a permutation
        ! vector to alter the default ordering in each case. The first face
        ! index is assigned to the node which is closest to the face corner A
        ! having the smallest global index. The next indices are created in 
        ! the direction of the face edge AB, with B the smallest possible 
        ! global index.
        !
        Perm = 0

        MinEntryInd = MINLOC(FaceNodes(1:3))
        SELECT CASE(MinEntryInd(1))
        CASE(1)
          IF (FaceNodes(3) < FaceNodes(2)) THEN
            Perm = (/ 1,3,2 /)
          ELSE
            Perm = (/ 1,2,3 /)
          END IF

        CASE(2)
          IF (FaceNodes(3) < FaceNodes(1)) THEN
            Perm = (/ 2,3,1 /)
          ELSE
            Perm = (/ 2,1,3 /)
          END IF          

        CASE(3)
          IF (FaceNodes(1) < FaceNodes(2)) THEN
            Perm = (/ 3,1,2 /)
          ELSE
            Perm = (/ 3,2,1 /)
          END IF
        END SELECT

      CASE DEFAULT
        CALL Fatal('LagrangeTriFacePermutation', &
            'For triangular faces Lagrange Element Degree < 5 supported currently')

      END SELECT
    END FUNCTION LagrangeTriFacePermutation



!------------------------------------------------------------------------------
  END FUNCTION GetLagrangeIndexes
!------------------------------------------------------------------------------   


 !> Find a representative DG index for a node index. Note that
 !> there may be several possibilities and this is just one of them.
 !------------------------------------------------------------------  
   FUNCTION NodeToDGIndex(Mesh,nodeind) RESULT ( dgind )

    TYPE(Mesh_t) :: Mesh
    INTEGER :: nodeind
    INTEGER :: dgind
    
    INTEGER :: i,j,t
    TYPE(Element_t), POINTER :: Element

    dgind = 0

    IF(nodeind < 1 ) THEN
      CALL Warn('NodeToDGIndex','Cannot find DG index for too small node index!')
      RETURN
    END IF
    IF(nodeind > Mesh % NumberOfNodes ) THEN
      CALL Warn('NodeToDGIndex','Cannot find DG index for too large node index!')
      RETURN
    END IF         
    
    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      DO i = 1,Element % TYPE % NumberOfNodes          
        IF( Element % NodeIndexes(i) == nodeind ) THEN
          IF(.NOT. ASSOCIATED( Element % DGIndexes ) ) THEN
            CALL Fatal('NodeToDGIndex','There are no DG indexes!')
          END IF
          dgind = Element % DGIndexes(i)
          EXIT
        END IF
      END DO
      IF(dgind > 0 ) EXIT
    END DO
    
  END FUNCTION NodeToDGIndex



!------------------------------------------------------------------------------
!> Split a mesh at zero levelset by adding new nodes at the interface.
!> The idea is to be able to better represent shapes that are not initially
!> presented by body fitted finite element mesh. 
!------------------------------------------------------------------------------
  FUNCTION SplitMeshLevelset(Mesh,Vlist) RESULT( NewMesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Vlist    
    TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: phi(:)
    INTEGER, ALLOCATABLE :: EdgeSplit(:)
    LOGICAL, ALLOCATABLE :: CutNode(:)
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: SplitReady    
    REAL(KIND=dp), POINTER :: u(:),v(:),w(:),x(:),y(:),z(:)
    REAL(KIND=dp) :: Eps
    INTEGER, POINTER :: NodeIndexes(:), EdgeIndexes(:)    
    INTEGER :: i, j, j2, j3, k, k2, k3, l, l2, l3, m, n, &
        n_old, n_new, n_cut, n_split, n_neg, n_pos
    INTEGER :: NoHits, NewElCnt, BCCnt, prevl, &
        NodeCnt, FaceCnt, Node, ParentId 
    LOGICAL :: Found, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Eptr,Parent 
    CHARACTER(LEN=MAX_NAME_LEN) :: str       
    INTEGER, POINTER :: Child(:,:)
    REAL(KIND=dp) :: h1,h2,hprod,r,s1,s2 
    REAL(KIND=dp), POINTER :: stime(:)
    INTEGER :: ierr, ParTmp(6), ParSizes(6)
    INTEGER :: BodyOffset, SgnNode, BodyCount, LevelsetBC
    LOGICAL :: PosOffset, BulkParent, Parallel
    CHARACTER(*), PARAMETER :: Caller = 'SplitMeshLevelset'
        
!------------------------------------------------------------------------------
    CALL Info( Caller, 'Splitting finite element mesh at zero levelset!', Level = 5 )

    IF ( .NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Warn(Caller,'Original mesh not associated!')
      RETURN
    END IF
        
    CALL ResetTimer(Caller)
    
    DO i=1,Mesh % NumberOfBulkElements
      n = Mesh % Elements(i) % TYPE % ElementCode 
      IF( n /= 303 .AND. n /= 504 ) THEN
        CALL Fatal(Caller,'Only linear triangles and tets can be split: '//TRIM(I2S(n)))
      END IF
    END DO
    
    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. Mesh % SingleMesh )
        
    CALL Info( Caller, '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( Caller, Message, Level=6 )

    ! At this stage the coordinates have not been added as variable.
    ! We cannot use the UDF's if these are not available. Also time
    ! is needed by default in some calls. 
    Var => VariableGet( Mesh % Variables,'time')
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      CALL VariableAdd( Mesh % Variables, Mesh, &
          Name='Coordinate 1',DOFs=1,Values=Mesh % Nodes % x )   
      CALL VariableAdd(Mesh % Variables,Mesh, &
          Name='Coordinate 2',DOFs=1,Values=Mesh % Nodes % y )    
      CALL VariableAdd(Mesh % Variables,Mesh, &
          Name='Coordinate 3',DOFs=1,Values=Mesh % Nodes % z )    
      ALLOCATE(stime(1)); stime(1) = 0.0_dp
      CALL VariableAdd( Mesh % Variables, Mesh, &
          Name='Time',DOFs=1, Values=sTime )
      CurrentModel % Variables => Mesh % Variables
    END IF
    
    ! Initialize the levelset function for all nodes
    n_old = Mesh % NumberOfNodes
    ALLOCATE( Phi(n_old) )
    
    str = ListGetString( Vlist,'Levelset Variable', Found)
    IF( Found ) THEN
      Var => VariableGet(Mesh % Variables, str)
      IF(.NOT. ASSOCIATED(Var) ) THEN
        CALL Fatal(Caller,'"Levelset Variable" requested, but not available: '//TRIM(str))
      END IF
      Phi = 1.0_dp
      ! We revert to nodal indexes since it will be easier in the future!
      DO i=1,n_old
        j = Var % Perm(i)
        IF(j>0) Phi(i) = Var % Values(j)
      END DO
    ELSE      
      DO i=1,n_old
        Phi(i) = ListGetRealAtNode(Vlist,'Levelset Function', i, Found)
        IF(.NOT. Found ) THEN
          CALL Fatal(Caller,'"Levelset Function" needed to enrich the mesh!')             
        END IF
      END DO
    END IF
    
    Eps = ListGetCReal( Vlist,'Levelset Epsilon',Found )
    IF(.NOT. Found ) Eps = 1.0e-3
    
    n_pos = COUNT( Phi > 0.0 )
    n_neg = COUNT( Phi < 0.0 ) 
        
    BodyOffset = ListGetInteger( Vlist,'Levelset Body Offset',Found ) 
    PosOffset = ListGetLogical( Vlist,'Levelset Offset Positive',Found ) 
    LevelsetBC = ListGetInteger( Vlist,'Levelset Boundary',Found )
    IF(.NOT. Found) LevelsetBC = CurrentModel % NumberOfBCs
    
    IF( Parallel ) THEN
      n_pos = ParallelReduction(n_pos) 
      n_neg = ParallelReduction(n_neg)
    END IF
    
    CALL Info(Caller,'Positive and negative values: '&
        //TRIM(I2S(n_pos))//' vs. '//TRIM(I2S(n_neg)),Level=7)    
    
    IF( n_pos == 0 .OR. n_neg == 0 ) THEN
      CALL Warn(Caller,'Nothing to do, no zero levelset available!')
      RETURN
    END IF
       
    ! We need edges in order to do the splitting!
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT. EdgesPresent) CALL FindMeshEdges( Mesh )
        
    ALLOCATE( EdgeSplit(Mesh % NumberOfEdges), CutNode(n_old) )    
    EdgeSplit = 0
    CutNode = .FALSE.
        
    j = 0
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      h1 = Phi(NodeIndexes(1))
      h2 = Phi(NodeIndexes(2))
      hprod = h1*h2
      IF( hprod < 0.0_dp ) THEN
        r = ABS(h2)/(ABS(h1)+ABS(h2))
        IF( r <= Eps ) THEN
          CutNode(NodeIndexes(2)) = .TRUE.
        ELSE IF(1.0-r < Eps ) THEN
          CutNode(NodeIndexes(1)) = .TRUE.
        ELSE
          j = j+1 
          EdgeSplit(i) = j
        END IF
      ELSE IF( ABS(hprod) < 1.0d-20 ) THEN
        IF(ABS(h1) < 1.0e-20) CutNode(NodeIndexes(1)) = .TRUE. 
        IF(ABS(h2) < 1.0e-20) CutNode(NodeIndexes(2)) = .TRUE.
      END IF
    END DO
    
    n_new = j
    CALL Info(Caller,'Number of additional nodes: '//TRIM(I2S(n_new)),Level=6)

    j = COUNT( CutNode )
    CALL Info(Caller,'Number of cut nodes: '//TRIM(I2S(j)),Level=6)
    
!   Update nodal coordinates:
!   -------------------------
    NodeCnt = n_old + n_new 

!   Create the new mesh
!   -------------------------------
    NewMesh => AllocateMesh()    
    NewMesh % SingleMesh = Mesh % SingleMesh
    NewMesh % Name = Mesh % Name   

    CALL AllocateVector( NewMesh % Nodes % x, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % y, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % z, NodeCnt )

!   shortcuts (u,v,w) old mesh  nodes,
!   (x,y,z) new mesh nodes:
!   ----------------------------------
    u => Mesh % Nodes % x
    v => Mesh % Nodes % y
    w => Mesh % Nodes % z

    x => NewMesh % Nodes % x
    y => NewMesh % Nodes % y
    z => NewMesh % Nodes % z
!
!   new mesh includes old mesh nodes:
!   ----------------------------------
    x(1:n_old) = u
    y(1:n_old) = v
    z(1:n_old) = w

!   add new nodes where edges are split:
!   ------------------------------------
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      j = EdgeSplit(i)
      IF( j > 0 ) THEN
        j = j + n_old
        h1 = Phi(NodeIndexes(1))
        h2 = Phi(NodeIndexes(2))
        r = ABS(h2)/(ABS(h1)+ABS(h2))
        x(j) = r*u(NodeIndexes(1)) + (1-r)*u(NodeIndexes(2))
        y(j) = r*v(NodeIndexes(1)) + (1-r)*v(NodeIndexes(2))
        z(j) = r*w(NodeIndexes(1)) + (1-r)*w(NodeIndexes(2))
      END IF
    END DO

    CALL Info(Caller,'Added new nodes on the splitted edges.', Level=10 )  

    
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MinEdgeDOFs = Mesh % MinEdgeDOFs
    NewMesh % MinFaceDOFs = Mesh % MinFaceDOFs
    NewMesh % MaxEdgeDOFs = Mesh % MaxEdgeDOFs
    NewMesh % MaxFaceDOFs = Mesh % MaxFaceDOFs
    NewMesh % MaxElementDOFs = Mesh % MaxElementDOFs
    NewMesh % MeshDim = Mesh % MeshDim

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt

!   Update bulk elements:
!   =====================
!
!   First count maximum number of new elements:
!   -------------------------------------------
    NewElCnt = 0
    BodyCount = 0
    DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Eold => Mesh % Elements(i)
      j = 1

      Found = .FALSE.
      IF( ASSOCIATED( Eold % EdgeIndexes ) ) THEN
        Found = ANY(EdgeSplit(Eold % EdgeIndexes) > 0 )
      ELSE
        CALL Fatal(Caller,'No edges for element: '//TRIM(I2S(i)))
      END IF
      
      IF( Found ) THEN
        SELECT CASE( Eold % TYPE % ElementCode/100 )                
        CASE(2)
          j = 2
        CASE(3)
          j = 3
        CASE(5)
          j = 6
        END SELECT
        ! There will also be additional BC elements on the cut!
        j = j + 1
      END IF
      NewElCnt = NewElCnt + j
    END DO
    
    CALL Info( Caller,'Maximum estimated count of new elements: '//TRIM(I2S(NewElCnt)), Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info(Caller,'New mesh allocated.', Level=20 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 6 )
    Child = 0
    CALL Info(Caller,'Array for bulk elements allocated.', Level=20 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes

!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)
       NodeIndexes => Eold % NodeIndexes       
       n = Eold % TYPE % NumberOfNodes                
       n_split = COUNT( EdgeSplit(Eold % EdgeIndexes) > 0 )

       ! We continue splitting until the element is exhausted
       SplitReady = .FALSE.

       ! Split elements to no more than 6 pieces
       DO m = 1,6
         NewElCnt = NewElCnt + 1
         Child(i,m) = NewElCnt
         Enew => NewMesh % Elements(NewElCnt)

         Enew = Eold
         Enew % TYPE => Eold % TYPE
         Enew % BodyId = Eold % BodyId
         Enew % PartIndex = Eold % PartIndex
         Enew % ElementIndex = NewElCnt
         Enew % NDOFs = Eold % NDOFs
         Enew % EdgeIndexes => NULL()
         Enew % FaceIndexes => NULL()
         Enew % BoundaryInfo => NULL()
         
         CALL AllocateVector( ENew % NodeIndexes, n)
         
         IF( n_split == 0 ) THEN
           Enew % NodeIndexes = NodeIndexes
           DO j=1,n
             IF(.NOT. CutNode(NodeIndexes(j)) ) THEN
               ! This is a representative node that is used to determine the sign of the
               ! new elements in order to decide whether to add offset for body or not. 
               SgnNode = j
               EXIT
             END IF
           END DO
           
           SplitReady = .TRUE.
         ELSE           
           n_cut = COUNT( CutNode(NodeIndexes) )
       
           IF ( Eold % TYPE % ElementCode == 303 ) THEN         
             ! Split triangle to four triangles split on one or two edges
             !-----------------------------------------------------------
             IF( n_split == 2 ) THEN
               DO j=1,3
                 IF( EdgeSplit( Eold % EdgeIndexes(j) ) == 0 ) EXIT
               END DO
               j2 = MODULO(j,3)+1
               j3 = MODULO(j+1,3)+1

               IF( m == 1 ) THEN
                 ! There are two ways to split the triangle.
                 ! Choose the one with shorter diameter.
                 s1 = (x(NodeIndexes(j)) - x(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2 + &
                     (y(NodeIndexes(j)) - y(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2 + &
                     (z(NodeIndexes(j)) - z(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2
                 s2 = (x(NodeIndexes(j2)) - x(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2 + &
                     (y(NodeIndexes(j2)) - y(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2 + &
                     (z(NodeIndexes(j2)) - z(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2
                 Enew % NodeIndexes(1) = NodeIndexes(j)
                 Enew % NodeIndexes(2) = NodeIndexes(j2)                 
                 IF( s1 < s2 ) THEN
                   Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 ELSE
                   Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                   
                 END IF
                 SgnNode = j
               ELSE IF(m==2) THEN
                 IF( s1 < s2 ) THEN
                   Enew % NodeIndexes(1) = NodeIndexes(j)
                   SgnNode = j
                 ELSE
                   Enew % NodeIndexes(1) = NodeIndexes(j2)                   
                   SgnNode = j2
                 END IF
                 Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                
               ELSE IF(m==3) THEN
                 Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))
                 Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 Enew % NodeIndexes(3) = NodeIndexes(j3)
                 SgnNode = j3
                 SplitReady = .TRUE.
               END IF

             ELSE IF( n_split == 1 ) THEN
               DO j=1,3
                 IF( EdgeSplit( Eold % EdgeIndexes(j) ) > 0 ) EXIT
               END DO
               j2 = MODULO(j,3)+1
               j3 = MODULO(j+1,3)+1

               ! One cut result to splitted elements only if the opposing node is cut through
               IF( .TRUE. .OR. CutNode(NodeIndexes(j3)) ) THEN
                 IF(m==1) THEN
                   Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
                   Enew % NodeIndexes(2) = NodeIndexes(j2)
                   Enew % NodeIndexes(3) = NodeIndexes(j3)
                   IF( CutNode(NodeIndexes(j3)) ) THEN
                     SgnNode = j2
                   ELSE
                     SgnNode = j3
                   END IF
                 ELSE IF(m==2) THEN
                   Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
                   Enew % NodeIndexes(2) = NodeIndexes(j3)
                   Enew % NodeIndexes(3) = NodeIndexes(j)
                   IF( CutNode(NodeIndexes(j3)) ) THEN
                     SgnNode = j
                   ELSE
                     SgnNode = j3
                   END IF 
                   SplitReady = .TRUE.
                 END IF
               ELSE
                 Enew % NodeIndexes = NodeIndexes
                 SgnNode = j3
                 SplitReady = .TRUE.
               END IF
             ELSE
               CALL Fatal(Caller,'Triangle can only deal with 1 and 2 splits!')
             END IF
           ELSE              
             CALL Fatal(Caller,'Element type '//TRIM(I2S(Eold % TYPE % ElementCode))//&
                 ' not supported by the levelset splitter.')
           END IF
         END IF

         ! Set offset for inside/outside elements of the zero levelset.
         ! The SgnNode is a representative node the sign of which tells whether we are inside
         ! or outside. 
         IF( PosOffset ) THEN
           IF( Phi(NodeIndexes(SgnNode)) > 0.0 )  THEN
             Enew % BodyId = Enew % BodyId + BodyOffset
             BodyCount = BodyCount + 1
           END IF
         ELSE
           IF( Phi(NodeIndexes(SgnNode)) < 0.0 )  THEN
             Enew % BodyId = Enew % BodyId + BodyOffset            
             BodyCount = BodyCount + 1
           END IF
         END IF
         IF( SplitReady ) EXIT
       END DO
     END DO
     
!   Update new mesh element counts:
!   -------------------------------
    NewMesh % NumberOfBulkElements = NewElCnt
    
    CALL Info(Caller,'Number of elements inside: '//TRIM(I2S(BodyCount)),Level=7)
    
   
!   Update boundary elements:
!   ---------------------------------------------------

    BCCnt = 0
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      IF( i == Mesh % NumberOfBulkElements + 1 ) THEN
        CALL Info(Caller,'Number of boundary elements from bulk cuts: '//TRIM(I2S(BCCnt)))           
        BCCnt = 0
      END IF
     
      Eold => Mesh % Elements(i)
      NodeIndexes => Eold % NodeIndexes             
      BulkParent = ( i <= Mesh % NumberOfBulkElements )
      n_split = COUNT( EdgeSplit(Eold % EdgeIndexes) > 0 )
      n_cut = COUNT( CutNode(NodeIndexes) )

      ! Elements created from bulk cuts require some splits or cuts.
      ! Existing boundary elements remain even without cuts.
      IF( BulkParent ) THEN
        IF( n_split + n_cut <= 1 ) CYCLE
      END IF
  
      SplitReady = .FALSE.
      
      ! Each existing boundary element may be cut to several pieces
      ! For triangles this is just max two!
      DO m=1,10          
        BCCnt = BCCnt + 1
        NewElCnt = NewElCnt + 1
        IF( NewElCnt > SIZE( NewMesh % Elements ) ) THEN
          CALL Fatal(Caller,'Too few elements allocated: '//TRIM(I2S(NewElCnt)))
        END IF
       
        Enew => NewMesh % Elements(NewElCnt)
        
        ALLOCATE(Enew % BoundaryInfo)         
        Enew % PartIndex = Eold % PartIndex
        Enew % ElementIndex = NewElCnt
        
        n = 2
        Enew % TYPE => GetElementType(202)
        CALL AllocateVector( ENew % NodeIndexes, n)
        Enew % NDOFs = n
        Enew % EdgeIndexes => NULL()
        Enew % FaceIndexes => NULL()
                
        IF( BulkParent ) THEN
          ! There are the new boundary elements that come from splitting the mesh
          ! at zero levelset. Give the boundary a new index. 
          Enew % BoundaryInfo % Constraint = LevelsetBC
          
          IF ( Eold % TYPE % ElementCode == 303 ) THEN         
            IF( n_split == 2 ) THEN
              DO j=1,3
                IF( EdgeSplit( Eold % EdgeIndexes(j) ) == 0 ) EXIT
              END DO
              j2 = MODULO(j,3)+1
              j3 = MODULO(j+1,3)+1
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
              Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                   
            ELSE IF( n_split == 1 .AND. n_cut == 1) THEN
              DO j=1,3
                IF( EdgeSplit( Eold % EdgeIndexes(j) ) > 0 ) EXIT
              END DO
              !j2 = MODULO(j,3)+1
              !j3 = MODULO(j+1,3)+1                        
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
              DO j2=1,3
                IF( CutNode(NodeIndexes(j2)) ) EXIT
              END DO
              Enew % NodeIndexes(2) = NodeIndexes(j2)
            ELSE IF( n_cut == 2) THEN
              DO j=1,3
                IF( .NOT. CutNode(NodeIndexes(j) ) ) EXIT
              END DO
              j2 = MODULO(j,3)+1
              j3 = MODULO(j+1,3)+1                        
              Enew % NodeIndexes(1) = NodeIndexes(j2)
              Enew % NodeIndexes(2) = NodeIndexes(j3)
            ELSE
              CALL Fatal(Caller,'Can only deal with 2 or 1+1 splits!')
            END IF
          ELSE              
            CALL Fatal(Caller,'Element type '//TRIM(I2S(Eold % TYPE % ElementCode))//&
                ' not supported by the levelset splitting.')
          END IF          
          SplitReady = .TRUE.
          
        ELSE
          ! Each existing boundary element may be cut to several pieces
          Enew % BoundaryInfo = Eold % BoundaryInfo
          
          IF( n_split == 0 ) THEN
            ! If no edge is split the element stays as is
            Enew % NodeIndexes = Eold % NodeIndexes
            SplitReady = .TRUE.
            
          ELSE IF( Eold % TYPE % ElementCode == 202 ) THEN
            IF(m==1) THEN
              Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
              Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(1))
            ELSE IF(m==2) THEN
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(1))
              Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
              SplitReady = .TRUE.
            END IF
          ELSE
            CALL Fatal(Caller,'Cannot do this element yet!')
          END IF
        END IF

         
        prevl = 0
        DO k=1,2
          ! Pointer to the found left/right bulk element
          Eptr => NULL()

          IF( BulkParent ) THEN
            ! If the boundary results from splitting existing elements then
            ! the parent is the existing bulk elements. 
            Parent => Mesh % Elements(i)
          ELSE
            ! If boundary results from existing boundary elements then the potential
            ! parents are the children of the old parents. 
            IF( k==1 ) THEN
              Parent => Eold % BoundaryInfo % Left
            ELSE            
              Parent => Eold % BoundaryInfo % Right
            END IF
            IF(.NOT. ASSOCIATED(Parent)) CYCLE
          END IF

          ! Find the correct parent among the splitted children of the
          ! initial bulk elements. There may be 1 or several children. 
          DO k2 = 1, 6            
            l = Child( Parent % ElementIndex, k2 )
            IF(l==0) CYCLE
            NoHits = 0
            
            IF( BulkParent ) THEN
              IF( k==2 .AND. l == prevl ) CYCLE
            END IF

            IF(l == 0 .OR. l > SIZE( NewMesh % Elements) ) THEN
            ! This is left for debugging...
#if 1
              PRINT *,'Size Elements:',l, SIZE(NewMesh % Elements)
              PRINT *,'Child:',m,n_split,n_cut,k2,BulkParent
              PRINT *,'Parent index:',i,Parent % ElementIndex
              PRINT *,'Parents children',Child( Parent % ElementIndex, :)
              PRINT *,'Parent indexes:',Parent % NodeIndexes
              PRINT *,'Enew:',Enew % NodeIndexes
              PRINT *,'Eold:',Eold % NodeIndexes
              PRINT *,'Eold edges:',Eold % EdgeIndexes
              PRINT *,'Eold edge node indexes:',Mesh % Edges(Eold % EdgeIndexes(1) ) % NodeIndexes
              DO l2=1,6
                IF(Child( Parent % ElementIndex, l2) == 0 ) EXIT
                PRINT *,'Parent:',l2,NewMesh % &
                    Elements(Child( Parent % ElementIndex,l2)) %  NodeIndexes
              END DO
              PRINT *,'old node indexes:',Mesh % Elements(Parent % ElementIndex) %  NodeIndexes
              PRINT *,'old edge indexes:',Mesh % Elements(Parent % ElementIndex) %  EdgeIndexes
              PRINT *,'old cut indexes:',CutNode(Mesh % Elements(Parent % ElementIndex) %  NodeIndexes)
              PRINT *,'old split indexes:',EdgeSplit(Mesh % Elements(Parent % ElementIndex) %  EdgeIndexes)
#endif
              EXIT
            END IF

            Eptr => NewMesh % Elements(l)

            DO l2 = 1,Enew % Type % NumberOfNodes 
              DO l3 = 1, Eptr % TYPE % NumberOfNodes
                IF( Enew % NodeIndexes(l2) == Eptr % NodeIndexes(l3) ) THEN
                  NoHits = NoHits + 1
                  EXIT
                END IF
              END DO
            END DO
            
            IF( NoHits == n ) EXIT
          END DO
          
          IF( NoHits == n ) THEN
            IF( k==1) THEN
              prevl = l
              Enew % BoundaryInfo % Left => Eptr
            ELSE
              Enew % BoundaryInfo % Right => Eptr
            END IF
          ELSE
            IF(k==1) CALL Warn(Caller,'Could not find even 1 parant!')
          END IF
            
        END DO
       
       ! When we have created all the new boundary elements resulting from splitting
       ! the master element then proceed to next element. 
       IF(SplitReady) EXIT
     END DO
   END DO


!   Update new mesh element counts:
!   -------------------------------
   CALL Info(Caller,'Number of total elements: '//TRIM(I2S(NewElCnt)),Level=7)
    
!   Update new mesh boundary element counts:
!   ----------------------------------------
   NewMesh % NumberOfBoundaryElements = NewElCnt - &
       NewMesh % NumberOfBulkElements
   NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs
   NewMesh % MaxElementNodes = Mesh % MaxElementNodes
   
    
   CALL Info( Caller, '******** New mesh ********', Level=6 )
   WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
   CALL Info( Caller, Message, Level=6 )
   WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
   CALL Info( Caller, Message, Level=6 )
   WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
   CALL Info( Caller, Message, Level=6 )


   ! Information of the new system size, also in parallel
   !----------------------------------------------------------------------
   ParTmp(1) = Mesh % NumberOfNodes
   ParTmp(2) = Mesh % NumberOfBulkElements
   ParTmp(3) = Mesh % NumberOfBoundaryElements
   ParTmp(4) = NewMesh % NumberOfNodes
   ParTmp(5) = NewMesh % NumberOfBulkElements
   ParTmp(6) = NewMesh % NumberOfBoundaryElements
   
   IF( .FALSE. .AND. Parallel ) THEN
     CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
     
     CALL Info(Caller,'Information on parallel mesh sizes',Level=8)
     CALL Info(Caller,'Initial mesh has '//TRIM(I2S(ParSizes(1)))//' nodes',Level=8)
     CALL Info(Caller,'Initial mesh has '//TRIM(I2S(ParSizes(2)))//' bulk elements',Level=8)
     CALL Info(Caller,'Initial mesh has '//TRIM(I2S(ParSizes(3)))//' boundary elements',Level=8)
     CALL Info(Caller,'New mesh has '//TRIM(I2S(ParSizes(4)))//' nodes',Level=5)
     CALL Info(Caller,'New mesh has '//TRIM(I2S(ParSizes(5)))//' bulk elements',Level=5)
     CALL Info(Caller,'New mesh has '//TRIM(I2S(ParSizes(6)))//' boundary elements',Level=5)
   END IF

    
   ! Update structures needed for parallel execution:
   !--------------------------------------------------
   IF( Parallel ) THEN
     CALL UpdateParallelInfo( Mesh, NewMesh )
   END IF

   ! Finalize:
   !-----------
   IF(.NOT.EdgesPresent) THEN
     CALL ReleaseMeshEdgeTables( Mesh )
     CALL ReleaseMeshFaceTables( Mesh )
   ELSE
     CALL FindMeshEdges( NewMesh )
   END IF

   CALL CheckTimer(Caller,Delete=.TRUE.)

   CALL Info(Caller,'Mesh was enriched with zero levelset',Level=8)

 CONTAINS
    
!------------------------------------------------------------------------------
    SUBROUTINE UpdateParallelInfo( Mesh, NewMesh )
!------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Edge
      INTEGER :: i,j1,j2,n,n0,m,istat
      LOGICAL :: Found
!------------------------------------------------------------------------------
!
!      Update mesh interfaces for parallel execution.
!      ==============================================
       n = NewMesh % NumberOfNodes
       ALLOCATE( NewMesh % ParallelInfo % NeighbourList(n), stat=istat )
       IF ( istat /= 0 ) CALL Fatal( Caller, 'Allocate error.' )
       DO i=1,n
         NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       CALL AllocateVector( NewMesh % ParallelInfo % NodeInterface,n  )       
       NewMesh % ParallelInfo % NodeInterface = .FALSE.

       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )
       NewMesh % ParallelInfo % GlobalDOFs = 0

       ! Inherit the old parallel data
       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % NodeInterface(1:n) = Mesh % ParallelInfo % NodeInterface
       NewMesh % ParallelInfo % GlobalDOFs(1:n) = Mesh % ParallelInfo % GlobalDOFs
       DO i=1,n
         m = SIZE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) 
         ALLOCATE( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours(m) )
         NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
             Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       END DO

       n0 = ParallelReduction(MAXVAL(Mesh % ParallelInfo % GlobalDofs),2)       
       CALL Info(Caller,'Offset for parallel numbering of new nodes: '//TRIM(I2S(n0)))

       ! We need global numbering for the edges that we use for the unique numbering of new nodes
       CALL SParEdgeNumbering(Mesh)
       
       DO i=1,Mesh % NumberOfEdges
         j = EdgeSplit(i)
         IF(j==0) CYCLE
         Edge => Mesh % Edges(j)

         ! Make a unique parallel number for the new nodes introduced at split edges
         NewMesh % ParallelInfo % GlobalDOFs(n+j) = n0 + Edge % GElementIndex         

         j1 = Edge % NodeIndexes(1)
         j2 = Edge % NodeIndexes(2)
         m = CountSameIntegers(Mesh % ParallelInfo % NeighbourList(j1) % Neighbours, &
             Mesh % ParallelInfo % NeighbourList(j2) % Neighbours, &
             NewMesh % ParallelInfo % NeighbourList(n+j) % Neighbours ) 
         NewMesh % ParallelInfo % NodeInterface(n+j) = (m>1)
       END DO
       
    END SUBROUTINE UpdateParallelInfo
    
  END FUNCTION SplitMeshLevelset
!------------------------------------------------------------------------------


  !> Create interface boundaries consisting of edges defined by the intersection of two higher
  !> dimensional boundaries. This may be useful for 3D meshes where 1D meshes have not been
  !> create in advance.
  !-------------------------------------------------------------------
  SUBROUTINE CreateIntersectionBCs(Model, Mesh)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element, Element2, Enew, Face, Face2, Parent
    INTEGER, POINTER :: NodeIndexes(:), NodeIndexes2(:), EdgeIndexes(:), EdgeIndexes2(:), ParentBCs(:)
    INTEGER :: i,i2,j,j2,k,k2,e,e2,l,n,n2,m,nbc,nbulk,nold,t,t2,istat,newbcs,newcnt,bc_id
    TYPE(Element_t), POINTER :: NewElements(:)
    TYPE(ValueList_t), POINTER :: BC
    INTEGER, ALLOCATABLE :: BoundaryId(:), IntersectionBCs(:,:)
    LOGICAL, ALLOCATABLE :: EdgeDone(:), NodeDone(:)
    LOGICAL :: Found, Hit, EdgesPresent 
    
    ! Count how many of the BCs are intersection BCs that we need to determine
    j = 0
    DO bc_id=1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values
      IF( ListCheckPresent( BC,'Intersection BC' ) ) j = j+1 
    END DO
    NewBCs = j
    IF(NewBCs==0) RETURN

    CALL Info('CreateIntersectionBCs',&
        'Number of intersection BCs to determine: '//TRIM(I2S(NewBCs)),Level=5)

    ! Create a fast look-up table that define the new BC indexes and the parent BCs
    ALLOCATE(IntersectionBCs(j,4))
    IntersectionBCs = 0
    j = 0
    DO bc_id=1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values
      ParentBCs => ListGetIntegerArray( BC,'Intersection BC',Found )
      IF(.NOT. Found ) CYCLE
      j = j + 1
      IF(SIZE(ParentBCs) /= 2 ) CALL Fatal('CreateIntersectionBCs','Only available for two parents!')
      IntersectionBCs(j,1) = Model % BCs(bc_id) % Tag
      IntersectionBCs(j,2:3) = ParentBCs(1:2)
    END DO

    nbulk = Mesh % NumberOfBulkElements
    nbc = Mesh % NumberOfBoundaryElements
    nold = nbulk + nbc

    ALLOCATE( BoundaryId( nbc ) )
    BoundaryId = 0

  
   
    DO t=1,nbc
      Element => Mesh % Elements(nbulk+t)
      
      ! Only treat 2D boundary elements
      IF( Mesh % MeshDim == 3 ) THEN
        IF( Element % TYPE % ElementCode < 300 ) CYCLE
      ELSE
        IF( Element % TYPE % ElementCode < 200 ) CYCLE
      END IF
        
      DO bc_id=1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
      END DO
      IF ( bc_id > Model % NumberOfBCs ) CYCLE

      IF( ANY(IntersectionBCs(:,2)==bc_id .OR. IntersectionBCs(:,3)==bc_id)) THEN
        BoundaryId(t) = bc_id
      END IF
    END DO

    n = COUNT( BoundaryId > 0 )
    CALL Info('CreateIntersectionBCs','Number of candidate intersection parents: '//TRIM(I2S(n)))
    
    ! Go the new boundary elements over two times.
    ! On the 1st loop just count the number of new elements.
    ! On the 2nd lopp add the new elements in the element list.
    !-------------------------------------------------------------    
    EdgesPresent = ASSOCIATED( Mesh % Edges )
    IF( Mesh % MeshDim == 3 ) THEN
      IF(.NOT. EdgesPresent ) THEN
        CALL FindMeshEdges( Mesh ) 
      END IF
      ALLOCATE( EdgeDone( Mesh % NumberOfEdges ) )       
      CALL CreateIntersection3D(.TRUE.,NewCnt)
      IF(NewCnt==0) THEN
        CALL Info('CreateIntersectionBCs','Could not find any additional interface elements!')        
        GOTO 1
      END IF
      CALL CreateIntersection3D(.FALSE.,NewCnt)
    ELSE
      ALLOCATE( NodeDone( Mesh % NumberOfNodes ) ) 
      CALL CreateIntersection2D(.TRUE.,NewCnt)
      IF(NewCnt==0) THEN
        CALL Info('CreateIntersectionBCs','Could not find any additional interface elements!')
        GOTO 1
      END IF
      CALL CreateIntersection2D(.FALSE.,NewCnt)
    END IF
               
    IF( InfoActive(10) ) THEN
      DO i=1,newbcs
        CALL Info('CreateIntersectionBCs','New boundary '//TRIM(I2S(IntersectionBCs(i,1)))//&
            ' with '//TRIM(I2S(IntersectionBCs(i,4)))//' elements')
      END DO
    END IF

1   IF(Mesh % MeshDim == 3 .AND. .NOT. EdgesPresent ) THEN
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    END IF
      
    CALL Info('CreateIntersectionBCs','All done!',Level=10)


  CONTAINS

    ! Find intersection between 2D boundaries i.e. the result will
    ! be a new 1D boundary.
    !-------------------------------------------------------------
    SUBROUTINE CreateIntersection3D(AllocateOnly,NewCnt)
      LOGICAL :: AllocateOnly
      INTEGER :: NewCnt
      
      EdgeDone = .FALSE.
      NewCnt = 0

      DO t=1,nbc
        j = BoundaryId(t) 
        IF(j==0) CYCLE
        Element => Mesh % Elements(nbulk+t)
        NodeIndexes => Element % NodeIndexes
        n = Element % Type % NumberOfNodes

        DO t2=t+1,nbc
          j2 = BoundaryId(t2) 
          IF(j2==0) CYCLE
          IF(j==j2) CYCLE
          Element2 => Mesh % Elements(nbulk+t2)
          NodeIndexes2 => Element2 % NodeIndexes
          n2 = Element2 % TYPE % NumberOfNodes

          ! Do we have any common nodes. Some are required...
          k = 0
          DO i=1,n
            IF( ANY(NodeIndexes(i) == NodeIndexes2(1:n2) ) ) k = k+1
          END DO
          IF(k<2) CYCLE
          
          DO l=1,newbcs
            ! Do we have a suitable pair of indexes for the parents
            IF( .NOT. ( ( IntersectionBCs(l,2) == j .AND. IntersectionBCs(l,3) == j2 ) .OR. &
                ( IntersectionBCs(l,2) == j2 .AND. IntersectionBCs(l,3) == j ) ) ) CYCLE
            
            EdgeIndexes => Element % EdgeIndexes
            IF(ASSOCIATED(EdgeIndexes)) THEN
              Face => Element
            ELSE
              Face => NULL()
              IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
                Face => Find_Face(Mesh, Element % BoundaryInfo % Left, Element )                
              END IF
              IF(.NOT. ASSOCIATED(Face) ) THEN
                IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
                  Face => Find_Face(Mesh, Element % BoundaryInfo % Right, Element )
                END IF
              END IF
              IF(ASSOCIATED( Face ) ) THEN
                EdgeIndexes => Face % EdgeIndexes
              ELSE
                CALL Fatal('CreateIntersectionBCs','EdgeIndexes not associated!')
              END IF
            END IF

            ! This is a probably candidate as we have two 2D elements of proper type
            ! sharing at least two nodes. Just have to find for which edges the intersection
            ! applies. It could sometimes be a false positive also. 
            DO i=1,Face % TYPE % NumberOfEdges
              e = EdgeIndexes(i)          
              IF( EdgeDone(e) ) CYCLE

              EdgeIndexes2 => Element2 % EdgeIndexes
              IF(ASSOCIATED(EdgeIndexes2)) THEN
                Face2 => Element2
              ELSE
                Face2 => NULL()
                IF( ASSOCIATED( Element2 % BoundaryInfo % Left ) ) THEN
                  Face2 => Find_Face(Mesh, Element2 % BoundaryInfo % Left, Element2 )
                END IF
                IF(.NOT. ASSOCIATED(Face2) ) THEN
                  IF( ASSOCIATED( Element2 % BoundaryInfo % Right ) ) THEN
                    Face2 => Find_Face(Mesh, Element2 % BoundaryInfo % Right, Element2 )
                  END IF
                END IF
                IF(ASSOCIATED( Face2 ) ) THEN
                  EdgeIndexes2 => Face2 % EdgeIndexes
                ELSE
                  CALL Fatal('CreateIntersectionBCs','EdgeIndexes2 not associated!')
                END IF
              END IF

              DO i2=1,Face2 % TYPE % NumberOfEdges
                e2 = EdgeIndexes2(i2)

                ! Ok, we have a hit. Same edge appearing in the proper parent
                ! boundary elements. Create the actual boundary element only if
                ! we have already allocated for it. 
                IF(e==e2) THEN
                  EdgeDone(e) = .TRUE.
                  NewCnt = NewCnt + 1

                  IF(.NOT. AllocateOnly ) THEN
                    Enew => Mesh % Elements(nold+NewCnt)        
                    ALLOCATE(Enew % BoundaryInfo)         
                    Enew % PartIndex = Element % PartIndex
                    Enew % ElementIndex = nold + NewCnt
                    Enew % TYPE => Mesh % Edges(e) % TYPE

                    m = Enew % TYPE % NumberOfNodes
                    CALL AllocateVector( ENew % NodeIndexes, m)
                    Enew % NodeIndexes = Mesh % Edges(e) % NodeIndexes
                    Enew % NDOFs = m
                    Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
                    Enew % BoundaryInfo % Left => Element
                    Enew % BoundaryInfo % Right => Element2

                    Enew % EdgeIndexes => NULL()
                    Enew % FaceIndexes => NULL()
                    Enew % PDefs => NULL()
                    Enew % BubbleIndexes => NULL()

                    ! Just a simple counter for the new BCs of this type
                    IntersectionBCs(l,4) = IntersectionBCs(l,4) + 1
                  END IF

                  EXIT
                END IF
              END DO

            END DO
          END DO
        END DO
      END DO

      ! There is nothing to do since no new elements will be created.
      IF( NewCnt == 0 ) RETURN
              
      IF(AllocateOnly) THEN
        ALLOCATE( NewElements(nold + NewCnt ) )
        CALL Info('CreateIntersectionBCs','Allocated for '//TRIM(I2S(NewCnt))//' new 1D boundary elements!',Level=6)
        
        NewElements(1:nold) = Mesh % Elements(1:nold)

        DO i=nbulk+1,nold
          Element => Mesh % Elements(i)        
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

          Parent => Element % BoundaryInfo % Left
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
          END IF

          Parent => Element % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
          END IF
        END DO

        DEALLOCATE(Mesh % Elements)
        Mesh % Elements => NewElements
        Mesh % NumberOfBoundaryElements = nbc + NewCnt
      END IF

    END SUBROUTINE CreateIntersection3D


    ! Find intersection between 1D boundaries i.e. the result will
    ! be a new 0D boundary (=node). 
    !-------------------------------------------------------------
    SUBROUTINE CreateIntersection2D(AllocateOnly,NewCnt)
      LOGICAL :: AllocateOnly
      INTEGER :: NewCnt
      
      NodeDone = .FALSE.
      NewCnt = 0

      DO t=1,nbc
        j = BoundaryId(t) 
        IF(j==0) CYCLE
        Element => Mesh % Elements(nbulk+t)
        NodeIndexes => Element % NodeIndexes
        n = Element % Type % NumberOfNodes

        DO t2=t+1,nbc
          j2 = BoundaryId(t2) 
          IF(j==j2 .OR. j2==0) CYCLE
          Element2 => Mesh % Elements(nbulk+t2)
          NodeIndexes2 => Element2 % NodeIndexes
          n2 = Element2 % TYPE % NumberOfNodes

          
          ! Do we have any common nodes. Some are required...
          k = 0
          e = 0
          DO i=1,n
            IF( ANY(NodeIndexes(i) == NodeIndexes2(1:n2) ) ) THEN
              e = NodeIndexes(i)
              k = k+1
            END IF
          END DO
          IF(k/=1) CYCLE


          DO l=1,newbcs

            ! Do we have a suitable pair of indexes for the parents

            IF( .NOT. ( ( IntersectionBCs(l,2) == j .AND. IntersectionBCs(l,3) == j2 ) .OR. &
                ( IntersectionBCs(l,2) == j2 .AND. IntersectionBCs(l,3) == j ) ) ) CYCLE

            
            NodeDone(e) = .TRUE.
            NewCnt = NewCnt + 1
            
            IF(.NOT. AllocateOnly ) THEN
              Enew => Mesh % Elements(nold+NewCnt)        
              ALLOCATE(Enew % BoundaryInfo)         
              Enew % PartIndex = Element % PartIndex
              Enew % ElementIndex = nold + NewCnt
              Enew % TYPE => GetElementType(101)
              
              CALL AllocateVector( ENew % NodeIndexes, 1)
              Enew % NodeIndexes = e
              Enew % NDOFs = 1
              Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
              Enew % BoundaryInfo % Left => Element
              Enew % BoundaryInfo % Right => Element2
              
              Enew % EdgeIndexes => NULL()
              Enew % FaceIndexes => NULL()
              Enew % PDefs => NULL()
              Enew % BubbleIndexes => NULL()

              IntersectionBCs(l,4) = IntersectionBCs(l,4) + 1
            END IF
            
            EXIT
          END DO
        END DO
      END DO
      
      IF( NewCnt == 0 ) RETURN
              
      IF(AllocateOnly) THEN
        ALLOCATE( NewElements(nold + NewCnt ) )
        CALL Info('CreateIntersectionBCs','Allocated for '//TRIM(I2S(NewCnt))//' new 0D boundary elements!',Level=6)
        
        NewElements(1:nold) = Mesh % Elements(1:nold)
        
        DO i=nbulk+1,nold
          Element => Mesh % Elements(i)        
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

          Parent => Element % BoundaryInfo % Left
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
          END IF

          Parent => Element % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
          END IF
        END DO

        DEALLOCATE(Mesh % Elements)
        Mesh % Elements => NewElements
        Mesh % NumberOfBoundaryElements = nbc + NewCnt
      END IF

    END SUBROUTINE CreateIntersection2D
    
    
  END SUBROUTINE CreateIntersectionBCs

  
!------------------------------------------------------------------------------
END MODULE MeshUtils
!------------------------------------------------------------------------------

!> \}

