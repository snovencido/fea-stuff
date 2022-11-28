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
!>  Routines for allocating and initialization of mesh structures.
!------------------------------------------------------------------------------

MODULE MeshBasics
  
  USE Types
  USE Messages
  USE ElementUtils
  USE ParallelUtils
  IMPLICIT NONE

CONTAINS


!------------------------------------------------------------------------------
!> Allocated one single element. 
!------------------------------------------------------------------------------
   FUNCTION AllocateElement() RESULT( Element )
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

     ALLOCATE( Element, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateElement', 'Unable to allocate a few bytes of memory?' )
     Element % BDOFs    =  0
     Element % NDOFs    =  0
     Element % BodyId   = -1
     Element % Splitted =  0
     Element % hK = 0
     Element % ElementIndex = 0
     Element % StabilizationMk = 0
     NULLIFY( Element % TYPE )
     NULLIFY( Element % PDefs )
     NULLIFY( Element % BubbleIndexes )
     NULLIFY( Element % DGIndexes )
     NULLIFY( Element % NodeIndexes )
     NULLIFY( Element % EdgeIndexes )
     NULLIFY( Element % FaceIndexes )
     NULLIFY( Element % BoundaryInfo )
!------------------------------------------------------------------------------
   END FUNCTION AllocateElement
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
   SUBROUTINE AllocatePDefinitions(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     TYPE(Element_t) :: Element

     ! Sanity check to avoid memory leaks
     IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
        ALLOCATE(Element % PDefs, STAT=istat)
        IF ( istat /= 0) CALL Fatal('AllocatePDefinitions','Unable to allocate memory')
     ELSE
       CALL Info('AllocatePDefinitions','P element definitions already allocated',Level=32)
     END IF

     ! Initialize fields
     Element % PDefs % P = 0 
     Element % PDefs % TetraType = 0
     Element % PDefs % isEdge = .FALSE.
     Element % PDefs % pyramidQuadEdge = .FALSE.
     Element % PDefs % localNumber = 0
     Element % PDefs % GaussPoints = 0
!------------------------------------------------------------------------------
   END SUBROUTINE AllocatePDefinitions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AllocateBoundaryInfo(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     TYPE(Element_t) :: Element

     ALLOCATE(Element % BoundaryInfo, STAT=istat)
     IF ( istat /= 0) CALL Fatal('AllocateBoundaryInfo','Unable to allocate memory')

     Element % BoundaryInfo % Left => NULL()
     Element % BoundaryInfo % Right => NULL()
     Element % BoundaryInfo % GebhardtFactors => NULL()
     Element % BoundaryInfo % Constraint =  0

!------------------------------------------------------------------------------
   END SUBROUTINE AllocateBoundaryInfo
!------------------------------------------------------------------------------

!> Allocate mesh structure and return handle to it.
!------------------------------------------------------------------------------
   FUNCTION AllocateMesh(NumberOfBulkElements, NumberOfBoundaryElements, &
       NumberOfNodes, InitParallel ) RESULT(Mesh)
!------------------------------------------------------------------------------
     INTEGER, OPTIONAL :: NumberOfBulkElements, NumberOfBoundaryElements, NumberOfNodes
     LOGICAL, OPTIONAL :: InitParallel
     TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
     INTEGER :: istat, i, n
     CHARACTER(*), PARAMETER :: Caller = 'AllocateMesh'
     
     ALLOCATE( Mesh, STAT=istat )
     IF ( istat /= 0 ) CALL Fatal( Caller, 'Unable to allocate a few bytes of memory?' )

!    Nothing computed on this mesh yet!
!    ----------------------------------
     Mesh % SavesDone    = 0
     Mesh % OutputActive = .FALSE.

     Mesh % AdaptiveDepth = 0
     Mesh % Changed   = .FALSE. !  TODO: Change this sometime
     Mesh % Stabilize = .FALSE.
     Mesh % MeshTag = 1

     Mesh % Variables => NULL()
     Mesh % Parent => NULL()
     Mesh % Child => NULL()
     Mesh % Next => NULL()
     Mesh % RootQuadrant => NULL()
     Mesh % Edges => NULL()
     Mesh % Faces => NULL()
     Mesh % Projector => NULL()
     Mesh % NumberOfEdges = 0
     Mesh % NumberOfFaces = 0

     Mesh % NumberOfBulkElements = 0
     Mesh % NumberOfBoundaryElements = 0
     Mesh % Elements => NULL()
     
     Mesh % DiscontMesh = .FALSE.
     Mesh % SingleMesh  = .FALSE.
     Mesh % InvPerm => NULL()

     Mesh % MinFaceDOFs = 1000
     Mesh % MinEdgeDOFs = 1000
     Mesh % MaxNDOFs = 0
     Mesh % MaxFaceDOFs = 0
     Mesh % MaxEdgeDOFs = 0
     Mesh % MaxBDOFs = 0
     Mesh % MaxElementDOFs  = 0
     Mesh % MaxElementNodes = 0

     Mesh % ViewFactors => NULL()

     ALLOCATE( Mesh % Nodes, STAT=istat )
     IF ( istat /= 0 ) CALL Fatal( Caller, 'Unable to allocate a few bytes of memory?' )
     
     NULLIFY( Mesh % Nodes % x )
     NULLIFY( Mesh % Nodes % y )
     NULLIFY( Mesh % Nodes % z )
     Mesh % Nodes % NumberOfNodes = 0
     Mesh % NumberOfNodes = 0
       
     Mesh % NodesOrig => Mesh % Nodes
     NULLIFY( Mesh % NodesMapped )

     Mesh % EntityWeightsComputed = .FALSE.
     Mesh % BCWeight => NULL()
     Mesh % BodyForceWeight => NULL()
     Mesh % BodyWeight => NULL()
     Mesh % MaterialWeight => NULL()
    
     Mesh % ParallelInfo % NumberOfIfDOFs =  0        
     NULLIFY( Mesh % ParallelInfo % GlobalDOFs )
     NULLIFY( Mesh % ParallelInfo % NodeInterface )
     NULLIFY( Mesh % ParallelInfo % NeighbourList )     

     i = 0
     IF( PRESENT( NumberOfBulkElements ) ) THEN       
       Mesh % NumberOfBulkElements = NumberOfBulkElements
       i = i + 1
     END IF
     
     IF( PRESENT( NumberOfBoundaryElements ) ) THEN
       Mesh % NumberOfBoundaryElements = NumberOfBoundaryElements
       i = i + 1
     END IF

     IF( PRESENT( NumberOfNodes ) ) THEN
       Mesh % NumberOfNodes = NumberOfNodes
       i = i + 1
     END IF
     
     IF( i > 0 ) THEN
       IF( i < 3 ) CALL Fatal(Caller,'Either give all or no optional parameters!')
       CALL InitializeMesh( Mesh, InitParallel )         
     END IF       
     
!------------------------------------------------------------------------------
   END FUNCTION AllocateMesh
!------------------------------------------------------------------------------


   ! Initialize mesh structures after the size information has been 
   ! retrieved.
   !----------------------------------------------------------------
   SUBROUTINE InitializeMesh(Mesh, InitParallel)     
     TYPE(Mesh_t), POINTER :: Mesh
     LOGICAL, OPTIONAL :: InitParallel
     
     INTEGER :: i,j,k,NoElems,istat
     TYPE(Element_t), POINTER :: Element
     CHARACTER(*), PARAMETER :: Caller = 'InitializeMesh'
     LOGICAL :: DoParallel
     
     IF( Mesh % NumberOfNodes == 0 ) THEN
       CALL Warn(Caller,'Mesh has zero nodes!')
       RETURN
     ELSE
       CALL Info(Caller,'Number of nodes in mesh: '&
           //TRIM(I2S(Mesh % NumberOfNodes)),Level=8)
     END IF

     CALL Info(Caller,'Number of bulk elements in mesh: '&
         //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=8)        

     CALL Info(Caller,'Number of boundary elements in mesh: '&
         //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=8)        

     Mesh % Nodes % NumberOfNodes = Mesh % NumberOfNodes          

     NoElems = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

     IF( NoElems == 0 ) THEN
       CALL Fatal('InitializeMesh','Mesh has zero elements!')
     END IF

     Mesh % MaxElementDOFs  = 0
     Mesh % MinEdgeDOFs     = 1000
     Mesh % MinFaceDOFs     = 1000
     Mesh % MaxEdgeDOFs     = 0
     Mesh % MaxFaceDOFs     = 0
     Mesh % MaxBDOFs        = 0

     Mesh % DisContMesh = .FALSE.
     Mesh % DisContPerm => NULL()
     Mesh % DisContNodes = 0

     CALL Info(Caller,'Initial number of max element nodes: '&
         //TRIM(I2S(Mesh % MaxElementNodes)),Level=10) 

     ! Allocate the elements
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Elements, NoElems, Caller )

     DO j=1,NoElems        
       Element => Mesh % Elements(j)        

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

     ! Allocate the nodes
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Nodes % x, Mesh % NumberOfNodes, Caller )
     CALL AllocateVector( Mesh % Nodes % y, Mesh % NumberOfNodes, Caller )
     CALL AllocateVector( Mesh % Nodes % z, Mesh % NumberOfNodes, Caller )
     
     IF( .NOT. PRESENT( InitParallel ) ) RETURN
     IF( .NOT. InitParallel ) RETURN
     
     CALL Info( Caller,'Allocating parallel info',Level=12)
     
     ALLOCATE(Mesh % ParallelInfo % GlobalDOFs(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     ALLOCATE(Mesh % ParallelInfo % NodeInterface(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     ALLOCATE(Mesh % ParallelInfo % NeighbourList(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     DO i=1,Mesh % NumberOfNodes
       NULLIFY(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
     END DO
     
   END SUBROUTINE InitializeMesh


!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMesh( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Projector_t), POINTER :: Projector
    TYPE(Projector_t), POINTER :: Projector1
    TYPE(Variable_t), POINTER  :: Var, Var1
    INTEGER :: i,j,k
    LOGICAL :: GotIt
    REAL(KIND=dp), POINTER :: ptr(:)
!------------------------------------------------------------------------------
 
!    Deallocate mesh variables:
!    --------------------------

    CALL Info('ReleaseMesh','Releasing mesh variables',Level=15)
    CALL ReleaseVariableList( Mesh % Variables )
    Mesh % Variables => NULL()

!    Deallocate mesh geometry (nodes,elements and edges):
!    ----------------------------------------------------
    IF ( ASSOCIATED( Mesh % Nodes ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh nodes',Level=15)
      IF ( ASSOCIATED( Mesh % Nodes % x ) ) DEALLOCATE( Mesh % Nodes % x )
      IF ( ASSOCIATED( Mesh % Nodes % y ) ) DEALLOCATE( Mesh % Nodes % y )
      IF ( ASSOCIATED( Mesh % Nodes % z ) ) DEALLOCATE( Mesh % Nodes % z )
      DEALLOCATE( Mesh % Nodes )

      IF ( ASSOCIATED( Mesh % ParallelInfo % GlobalDOFs ) ) &
          DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs )

      IF ( ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) THEN 
        DO i=1,Mesh % NumberOfNodes
          IF(ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) &
              DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
        END DO
        DEALLOCATE( Mesh % ParallelInfo % NeighbourList )
      END IF

      IF ( ASSOCIATED( Mesh % ParallelInfo % NodeInterface ) ) &
          DEALLOCATE( Mesh % ParallelInfo % NodeInterface )
    END IF

    Mesh % Nodes => NULL()

    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh edges',Level=15)
      CALL ReleaseMeshEdgeTables( Mesh )
      Mesh % Edges => NULL()
    END IF

    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh faces',Level=15)
      CALL ReleaseMeshFaceTables( Mesh )
      Mesh % Faces => NULL()
    END IF

    IF (ASSOCIATED(Mesh % ViewFactors) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh view factors',Level=15)
      CALL ReleaseMeshFactorTables( Mesh % ViewFactors )
      Mesh % ViewFactors => NULL()
    END IF


!    Deallocate mesh to mesh projector structures:
!    ---------------------------------------------
    Projector => Mesh % Projector
    DO WHILE( ASSOCIATED( Projector ) )
      CALL Info('ReleaseMesh','Releasing mesh projector',Level=15)
      CALL FreeMatrix( Projector % Matrix )
      CALL FreeMatrix( Projector % TMatrix )
      Projector1 => Projector
      Projector => Projector % Next
      DEALLOCATE( Projector1 )
    END DO
    Mesh % Projector => NULL()


!    Deallocate quadrant tree (used in mesh to mesh interpolation):
!    --------------------------------------------------------------
    IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh quadrant tree',Level=15)
      CALL FreeQuadrantTree( Mesh % RootQuadrant )
      Mesh % RootQuadrant => NULL()
    END IF


    IF ( ASSOCIATED( Mesh % Elements ) ) THEN
      CALL Info('ReleaseMesh','Releasing mesh elements',Level=15)

      DO i=1,SIZE(Mesh % Elements)

!          Boundaryinfo structure for boundary elements
!          ---------------------------------------------
        IF ( Mesh % Elements(i) % Copy ) CYCLE

        IF ( i > Mesh % NumberOfBulkElements ) THEN
          IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo ) ) THEN
            IF (ASSOCIATED(Mesh % Elements(i) % BoundaryInfo % GebhardtFactors)) THEN
              IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo % &
                  GebhardtFactors % Elements ) ) THEN
                DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                    GebhardtFactors % Elements )
                DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                    GebhardtFactors % Factors )
              END IF
              DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % GebhardtFactors )
            END IF
            DEALLOCATE( Mesh % Elements(i) % BoundaryInfo )
          END IF
        END IF

        IF ( ASSOCIATED( Mesh % Elements(i) % NodeIndexes ) ) &
            DEALLOCATE( Mesh % Elements(i) % NodeIndexes )
        Mesh % Elements(i) % NodeIndexes => NULL()

        IF ( ASSOCIATED( Mesh % Elements(i) % DGIndexes ) ) &
            DEALLOCATE( Mesh % Elements(i) % DGIndexes )
        Mesh % Elements(i) % DGIndexes => NULL()

        IF ( ASSOCIATED( Mesh % Elements(i) % BubbleIndexes ) ) &
            DEALLOCATE( Mesh % Elements(i) % BubbleIndexes )
        Mesh % Elements(i) % BubbleIndexes => NULL()

        ! This creates problems later on!!!
        !IF ( ASSOCIATED( Mesh % Elements(i) % PDefs ) ) &
        !   DEALLOCATE( Mesh % Elements(i) % PDefs )
        
        Mesh % Elements(i) % PDefs => NULL() 
      END DO
      
      DEALLOCATE( Mesh % Elements )
      Mesh % Elements => NULL()
    END IF
     
    Mesh % NumberOfNodes = 0
    Mesh % NumberOfBulkElements = 0
    Mesh % NumberOfBoundaryElements = 0
    
    CALL Info('ReleaseMesh','Releasing mesh finished',Level=15)
    
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshEdgeTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
      CALL Info('ReleaseMeshEdgeTables','Releasing number of edges: '&
          //TRIM(I2S(Mesh % NumberOfEdges)),Level=30)
      
       DO i=1,Mesh % NumberOfEdges
          Edge => Mesh % Edges(i)
          IF ( ASSOCIATED( Edge % NodeIndexes ) ) THEN
             DEALLOCATE( Edge % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Edge % BoundaryInfo ) ) THEN
             DEALLOCATE( Edge % BoundaryInfo )
          END IF
       END DO
       DEALLOCATE( Mesh % Edges )

       NULLIFY( Mesh % Edges )
       IF( Mesh % NumberOfEdges == 0 ) RETURN
       Mesh % NumberOfEdges = 0
       
       IF( ASSOCIATED( Mesh % Elements ) ) THEN      
         DO i=1,SIZE(Mesh % Elements)
           IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) THEN
             DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
             Mesh % Elements(i) % EdgeIndexes => NULL()
           END IF
         END DO
       END IF
     END IF
       
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshEdgeTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFaceTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Face
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
      CALL Info('ReleaseMeshFaceTables','Releasing number of faces: '&
          //TRIM(I2S(Mesh % NumberOfFaces)))

      DO i=1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i)
          IF ( ASSOCIATED( Face % NodeIndexes ) ) THEN
             DEALLOCATE( Face % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Face % BoundaryInfo ) ) THEN
             DEALLOCATE( Face % BoundaryInfo )
          END IF
       END DO

       DEALLOCATE( Mesh % Faces )
       NULLIFY( Mesh % Faces )
       IF( Mesh % NumberOfFaces == 0 ) RETURN
       
       Mesh % NumberOfFaces = 0

       IF( ASSOCIATED( Mesh % Elements ) ) THEN
         DO i=1,SIZE(Mesh % Elements)
           IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) THEN
             DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
             Mesh % Elements(i) % FaceIndexes => NULL()
           END IF
         END DO
       END IF
     END IF
       
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFaceTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFactorTables( Factors )
!------------------------------------------------------------------------------
    TYPE(Factors_t), POINTER :: Factors(:)
!------------------------------------------------------------------------------
    INTEGER :: i
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Factors ) ) THEN
       DO i=1,SIZE( Factors)
          IF (ASSOCIATED(Factors(i) % Factors))  DEALLOCATE(Factors(i) % Factors)
          IF (ASSOCIATED(Factors(i) % Elements)) DEALLOCATE(Factors(i) % Elements)
       END DO
       DEALLOCATE(  Factors )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFactorTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION Find_Edge(Mesh,Parent,Element) RESULT(ptr)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Ptr
    TYPE(Mesh_t) :: Mesh
    TYPE(Element_t) :: Parent, Element

    INTEGER :: i,j,k,n

    Ptr => NULL()
    DO i=1,Parent % TYPE % NumberOfEdges
      Ptr => Mesh % Edges(Parent % EdgeIndexes(i))
      n=0
      DO j=1,Ptr % TYPE % NumberOfNodes
        DO k=1,Element % TYPE % NumberOfNodes
          IF (Ptr % NodeIndexes(j) == Element % NodeIndexes(k)) n=n+1
        END DO
      END DO
      IF (n==Ptr % TYPE % NumberOfNodes) EXIT
    END DO
!------------------------------------------------------------------------------
  END FUNCTION Find_Edge
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
  FUNCTION Find_Face(Mesh,Parent,Element) RESULT(ptr)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Ptr
    TYPE(Mesh_t) :: Mesh
    TYPE(Element_t) :: Parent, Element

    INTEGER :: i,j,k,n

    Ptr => NULL()
    DO i=1,Parent % TYPE % NumberOfFaces
      Ptr => Mesh % Faces(Parent % FaceIndexes(i))
      n=0
      DO j=1,Ptr % TYPE % NumberOfNodes
        DO k=1,Element % TYPE % NumberOfNodes
          IF (Ptr % NodeIndexes(j) == Element % NodeIndexes(k)) n=n+1
        END DO
      END DO
      IF (n==Ptr % TYPE % NumberOfNodes) EXIT
    END DO
!------------------------------------------------------------------------------
  END FUNCTION Find_Face
!------------------------------------------------------------------------------


  !> Returns the local nodal coordinate values from the global mesh
  !> structure in the given Element and Indexes.
  !---------------------------------------------------------------------------
  SUBROUTINE CopyElementNodesFromMesh( ElementNodes, Mesh, n, Indexes)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Mesh_t) :: Mesh
    INTEGER :: n,m
    INTEGER, POINTER :: Indexes(:)

    IF ( .NOT. ASSOCIATED( ElementNodes % x ) ) THEN
      m = n
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n),ElementNodes % z(n) )
    ELSE
      m = SIZE(ElementNodes % x)
      IF ( m < n ) THEN
        DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z)
        ALLOCATE( ElementNodes % x(n), ElementNodes % y(n),ElementNodes % z(n) )
      ELSE IF( m > n ) THEN
        ElementNodes % x(n+1:m) = 0.0_dp
        ElementNodes % y(n+1:m) = 0.0_dp
        ElementNodes % z(n+1:m) = 0.0_dp
      END IF
    END IF

    ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
    ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
    ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))

  END SUBROUTINE CopyElementNodesFromMesh


  ! Mark nodes that are associated with at least some boundary element.
  !------------------------------------------------------------------------------
  SUBROUTINE MarkBCNodes(Mesh,BCNode,NoBCNodes)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, ALLOCATABLE :: BCNode(:)
    INTEGER :: NoBCNodes

    INTEGER :: elem
    TYPE(Element_t), POINTER :: Element

    CALL Info('MarkInterfaceNodes','Marking interface nodes',Level=8)

    IF(.NOT. ALLOCATED( BCNode ) ) THEN
      ALLOCATE( BCNode( Mesh % NumberOfNodes ) )
    END IF
    BCNode = .FALSE. 

    DO elem=Mesh % NumberOfBulkElements + 1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      Element => Mesh % Elements( elem )         
      !IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE

      BCNode(Element % NodeIndexes) = .TRUE.
    END DO

    NoBCNodes = COUNT( BCNode )

    CALL Info('MarkBCNodes','Number of BC nodes: '//TRIM(I2S(NoBCNodes)),Level=8)

  END SUBROUTINE MarkBCNodes

  
  
  !---------------------------------------------------------------------------
  ! Simply fitting of cylinder into a point cloud. This is done in two phases.
  ! 1) The axis of the cylinder is found by minimizing the \sum((n_i*t)^2)
  !    for each component of of t where n_i:s are the surface normals. 
  !    This is fully generic and assumes no positions. 
  ! 2) The radius and center point of the cylinder are found by fitting a circle
  !    in the chosen plane to three representative points. Currently the fitting
  !    can only be done in x-y plane. 
  !---------------------------------------------------------------------------
  SUBROUTINE CylinderFit(PMesh, PParams, BCind, dim, FitParams) 
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Valuelist_t), POINTER :: PParams
    INTEGER, OPTIONAL :: BCind
    INTEGER, OPTIONAL :: dim
    REAL(KIND=dp), OPTIONAL :: FitParams(:)
    
    INTEGER :: i,j,k,n,t,AxisI,iter,cdim,ierr
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: NiNj(9),A(3,3),F(3),M11,M12,M13,M14
    REAL(KIND=dp) :: d1,d2,MinDist,MaxDist,Dist,X0,Y0,Rad
    REAL(KIND=dp) :: Normal(3), AxisNormal(3), Tangent1(3), Tangent2(3), Coord(3), &
        CircleCoord(9)
    INTEGER :: CircleInd(3) 
    LOGICAL :: BCMode, DoIt, GotNormal, GotCenter, GotRadius
    INTEGER :: Tag, t1, t2
    LOGICAL, ALLOCATABLE :: ActiveNode(:)
    REAL(KIND=dp), POINTER :: rArray(:,:)
    

    BCMode = PRESENT( BCind )

    ! Set the range for the possible active elements. 
    ! Set the range for the possible active elements. 
    IF( BCMode ) THEN
      t1 = PMesh % NumberOfBulkElements + 1
      t2 = PMesh % NumberOfBulkElements + PMesh % NumberOfBoundaryElements
      Tag = CurrentModel % BCs(BCind) % Tag
      ALLOCATE( ActiveNode( PMesh % NumberOfNodes ) )
      ActiveNode = .FALSE.
    ELSE
      t1 = 1
      t2 = PMesh % NumberOfBulkElements      
    END IF
    
    ! If this is a line mesh there is really no need to figure out the 
    ! direction of the rotational axis. It can only be aligned with the z-axis.
    DO t=t1, t2
      Element => PMesh % Elements(t)
      IF( BCMode ) THEN
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE     
        IF ( Element % BoundaryInfo % Constraint /= Tag ) CYCLE
      END IF
      IF( Element % TYPE % ElementCode < 300 ) THEN
        cdim = 2
      ELSE
        cdim = 3
      END IF
      EXIT
    END DO
    
    IF( BcMode ) THEN
      cdim = ParallelReduction( cdim, 2 )
    END IF

    AxisNormal = 0.0_dp
    IF( cdim == 2 ) THEN
      GotNormal = .TRUE.
      AxisNormal(3) = 1.0_dp
    ELSE      
      rArray => ListGetConstRealArray( PParams,'Cylinder Normal',GotNormal)
      IF( GotNormal) AxisNormal(1:3) = rArray(1:3,1)
    END IF

    Coord = 0.0_dp
    rArray => ListGetConstRealArray( PParams,'Cylinder Center',GotCenter)
    IF( GotCenter) Coord(1:cdim) = rArray(1:cdim,1)
    
    Rad = ListGetConstReal( PParams,'Cylinder Radius',GotRadius)
 
    ! Do we have the fitting done already? 
    IF( GotNormal .AND. GotCenter .AND. GotRadius ) THEN
      IF( PRESENT(FitParams) ) THEN
        CALL Info('CylinderFit','Using cylinder paramaters from list',Level=25)
        FitParams(1:cdim) = Coord(1:cdim)
        IF( cdim == 2 ) THEN
          FitParams(3) = Rad
        ELSE
          FitParams(4:6) = AxisNormal
          FitParams(7) = Rad
        END IF
      END IF
      RETURN
    END IF
                  
    n = PMesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

       
    ! Compute the inner product of <N*N> for the elements
    NiNj = 0.0_dp
    DO t=t1, t2
      Element => PMesh % Elements(t)

      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes

      IF( BCMode ) THEN
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE     
        IF ( Element % BoundaryInfo % Constraint /= Tag ) CYCLE
        ActiveNode(Element % NodeIndexes(1:n)) = .TRUE.
      END IF
              
      ! If we know the Normal we only tag the boundary nodes
      IF(GotNormal) CYCLE

      Nodes % x(1:n) = PMesh % Nodes % x(NodeIndexes(1:n))
      Nodes % y(1:n) = PMesh % Nodes % y(NodeIndexes(1:n))
      Nodes % z(1:n) = PMesh % Nodes % z(NodeIndexes(1:n))           
      
      Normal = NormalVector( Element, Nodes, Check = .FALSE. ) 
      DO i=1,3
        DO j=1,3
          NiNj(3*(i-1)+j) = NiNj(3*(i-1)+j) + Normal(i) * Normal(j)
        END DO
      END DO
    END DO
      
    IF(GotNormal) GOTO 100 

    ! Only in BC mode we do currently parallel reduction.
    ! This could be altered too. 
    IF( BCMode ) THEN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,NiNj,9, &
          MPI_DOUBLE_PRECISION,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF
      
    ! The potential direction for the cylinder axis is the direction with 
    ! least hits for the normal.
    AxisI = 1 
    DO i=2,3
      IF( NiNj(3*(i-1)+i) < NiNj(3*(AxisI-1)+AxisI) ) AxisI = i 
    END DO

    CALL Info('CylinderFit','Axis coordinate set to be: '//TRIM(I2S(AxisI)))

    ! Keep the dominating direction fixed and iteratively solve the two other directions
    AxisNormal = 0.0_dp
    AxisNormal(AxisI) = 1.0_dp

    ! Basically we could solve from equation Ax=0 the tangent but only up to a constant.
    ! Thus we enforce the axis direction to one by manipulation the matrix equation 
    ! thereby can get a unique solution. 
    DO i=1,3
      DO j=1,3
        A(i,j) = NiNj(3*(i-1)+j)
      END DO
    END DO
    A(AxisI,1:3) = 0.0_dp
    A(AxisI,AxisI) = 1.0_dp
    CALL InvertMatrix( A, 3 )
    AxisNormal = A(1:3,AxisI)

    ! Normalize the axis normal length to one    
    AxisNormal = AxisNormal / SQRT( SUM( AxisNormal ** 2 ) )
    IF( 1.0_dp - MAXVAL( ABS( AxisNormal ) ) > 1.0d-5 ) THEN
      CALL Warn('CylinderFit','The cylinder axis is not aligned with any axis!')
    END IF

100 CALL TangentDirections( AxisNormal,Tangent1,Tangent2 )

    IF( InfoActive(30) .AND. ParEnv % MyPe == 0 ) THEN
      PRINT *,'Axis Normal:',AxisNormal
      PRINT *,'Axis Tangent 1:',Tangent1
      PRINT *,'Axis Tangent 2:',Tangent2
    END IF

    ! Finding three points with maximum distance in the tangent directions

    ! First, find the single extremum point in the first tangent direction
    ! Save the local coordinates in the N-T system of the cylinder
    MinDist = HUGE(MinDist)
    MaxDist = -HUGE(MaxDist)
    DO i=1, PMesh % NumberOfNodes
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)

      IF( BCMode ) THEN
        IF( .NOT. ActiveNode(i) ) CYCLE
      END IF
      
      d1 = SUM( Tangent1 * Coord )
      IF( d1 < MinDist ) THEN
        MinDist = d1
        CircleInd(1) = i
      END IF
      IF( d1 > MaxDist ) THEN
        MaxDist = d1
        CircleInd(2) = i
      END IF
    END DO

    CircleCoord = -HUGE(CircleCoord)
    DO j=1,2    
      i = CircleInd(j)

      IF( BCMode ) THEN
        IF(j==1) THEN
          Dist = ParallelReduction( MinDist, 1 )
          IF(ABS(MinDist-Dist) > 1.0e-10) CYCLE
        ELSE IF(j==2) THEN
          Dist = ParallelReduction( MaxDist, 2)
          IF(ABS(MaxDist-Dist) > 1.0e-10) CYCLE
        END IF
      END IF
        
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      CircleCoord(3*(j-1)+1) = SUM( Tangent1 * Coord ) 
      CircleCoord(3*(j-1)+2) = SUM( Tangent2 * Coord ) 
      CircleCoord(3*(j-1)+3) = SUM( AxisNormal * Coord )
    END DO

    IF( BCMode ) THEN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,CircleCoord,6, &
          MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
      
    !PRINT *,'MinDist1:',MinDist,CircleInd(1),CircleCoord(1,:)

    ! Find one more point such that their minimum distance to the previous point(s)
    ! is maximized. This takes some time but the further the nodes are apart the more 
    ! accurate it will be to fit the circle to the points. Also if there is just 
    ! a symmetric section of the cylinder it is important to find the points rigorously.
    j = 3
    ! The maximum minimum distance of any node from the previously defined nodes
    MaxDist = 0.0_dp
    DO i=1, PMesh % NumberOfNodes
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      IF( BCMode ) THEN
        IF( .NOT. ActiveNode(i) ) CYCLE
      END IF
      
      ! Minimum distance from the previously defined nodes
      MinDist = HUGE(MinDist)
      DO k=1,j-1
        d1 = SUM( Tangent1 * Coord )
        d2 = SUM( Tangent2 * Coord )
        Dist = ( d1 - CircleCoord(3*(k-1)+1) )**2 + ( d2 - CircleCoord(3*(k-1)+2) )**2
        MinDist = MIN( Dist, MinDist )
      END DO
      
      ! If the minimum distance to either previous selelected nodes
      ! is greater than in any other node, choose this
      IF( MaxDist < MinDist ) THEN
        MaxDist = MinDist 
        CircleInd(j) = i
      END IF
    END DO
    
    ! Ok, we have found the point now set the circle coordinates 
    DoIt = .TRUE.
    IF( BCMode ) THEN
      Dist = ParallelReduction( MaxDist, 2 )
      DoIt = ( ABS(MaxDist-Dist) < 1.0e-10 )
    END IF

    IF( DoIt ) THEN
      i = CircleInd(j)
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      CircleCoord(3*(j-1)+1) = SUM( Tangent1 * Coord ) 
      CircleCoord(3*(j-1)+2) = SUM( Tangent2 * Coord ) 
      CircleCoord(3*(j-1)+3) = SUM( AxisNormal * Coord )
    END IF

    IF( BCMode ) THEN
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,CircleCoord,9, &
          MPI_DOUBLE_PRECISION,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
      
    IF( InfoActive(30) .AND. ParEnv % MyPe == 0 ) THEN
      DO i=1,3
        PRINT *,'Circle Coord:',i,CircleCoord(3*i-2:3*i) 
      END DO
    END IF
      
    ! Given three nodes it is possible to analytically compute the center point and
    ! radius of the cylinder from a 4x4 determinant equation. The matrices values
    ! m1i are the determinants of the comatrices. 

    A(1:3,1) = CircleCoord(1::3)  ! x
    A(1:3,2) = CircleCoord(2::3)  ! y
    A(1:3,3) = 1.0_dp
    m11 = Det3x3( a )

    A(1:3,1) = CircleCoord(1::3)**2 + CircleCoord(2::3)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(2::3)  ! y
    A(1:3,3) = 1.0_dp
    m12 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1::3)**2 + CircleCoord(2::3)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(1::3)  ! x
    A(1:3,3) = 1.0_dp
    m13 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1::3)**2 + CircleCoord(2::3)**2 ! x^2+y^2
    A(1:3,2) = CircleCoord(1::3)  ! x
    A(1:3,3) = CircleCoord(2::3)  ! y
    m14 = Det3x3( a )

    !PRINT *,'determinants:',m11,m12,m13,m14

    IF( ABS( m11 ) < EPSILON( m11 ) ) THEN
      CALL Fatal('CylinderFit','Points cannot be an a circle')
    END IF

    X0 =  0.5_dp * m12 / m11 
    Y0 = -0.5_dp * m13 / m11
    rad = SQRT( x0**2 + y0**2 + m14/m11 )

    Coord = x0 * Tangent1 + y0 * Tangent2

    IF( InfoActive(30) .AND. ParEnv % MyPe == 0) THEN
      PRINT *,'Cylinder center and radius:',Coord, rad
    END IF

    ALLOCATE( rArray(3,1) )
    rArray(1:3,1) = Coord 
    CALL ListAddConstRealArray( PParams,'Cylinder Center', 3, 1, rArray ) 
    IF(.NOT. GotNormal ) THEN
      rArray(1:3,1) = AxisNormal 
      CALL ListAddConstRealArray( PParams,'Cylinder Normal', 3, 1, rArray ) 
    END IF
    DEALLOCATE( rArray ) 
    CALL ListAddConstReal( PParams,'Cylinder Radius',rad )

    IF( PRESENT( FitParams ) ) THEN
      IF( cdim == 2 ) THEN
        FitParams(1:2) = Coord(1:2)
        FitParams(3) = rad
      ELSE
        FitParams(1:3) = Coord(1:3)
        FitParams(4:6) = AxisNormal(1:3)
        FitParams(7) = rad
      END IF
    END IF
      
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

  CONTAINS
    
    ! Compute the value of 3x3 determinant
    !-------------------------------------------
    FUNCTION Det3x3( A ) RESULT ( val ) 
      
      REAL(KIND=dp) :: A(:,:)
      REAL(KIND=dp) :: val

      val = A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) &
          - A(1,2) * ( A(2,1) * A(3,3) - A(2,3) * A(3,1) ) &
          + A(1,3) * ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) 

    END FUNCTION Det3x3

  END SUBROUTINE CylinderFit



  ! Code for fitting a sphere. Not yet used.
  !-------------------------------------------------------------------------
  SUBROUTINE SphereFit(Mesh, Params, BCind, FitParams ) 
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, OPTIONAL :: BCind
    REAL(KIND=dp), OPTIONAL :: FitParams(:)

    INTEGER :: i,j,t,t1,t2,NoNodes,Tag
    LOGICAL :: BCMode
    LOGICAL, ALLOCATABLE :: ActiveNode(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)    
    REAL(KIND=dp) :: xc,yc,zc,Rad

    IF( PRESENT( FitParams ) ) THEN
      IF( ListCheckPresent( Params,'Sphere Radius') ) THEN
        CALL Info('SphereFit','Using predefined values for sphere parameters',Level=20)
        FitParams(1) = ListGetConstReal( Params,'Sphere Center X')
        FitParams(2) = ListGetConstReal( Params,'Sphere Center Y')
        FitParams(3) = ListGetConstReal( Params,'Sphere Center Z')
        FitParams(4) = ListGetConstReal( Params,'Sphere Radius')
        RETURN
      END IF
    END IF
      
    
    CALL Info('SphereFit','Trying to fit a sphere to element patch',Level=6)

    ! Set the range for the possible active elements. 
    IF( PRESENT( BCind ) ) THEN
      BCMode = .TRUE.
      t1 = Mesh % NumberOfBulkElements + 1
      t2 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Tag = CurrentModel % BCs(BCind) % Tag
      ALLOCATE( ActiveNode( Mesh % NumberOfNodes ) )
      ActiveNode = .FALSE.
    ELSE
      BCMode = .FALSE.
      t1 = 1
      t2 = Mesh % NumberOfBulkElements
    END IF

    ! Mark the nodes that belong to the active elements.
    ! 1) Either we only have bulk elements in which case we use all of the nodes or
    ! 2) We are given a boundary index and only use the nodes related to it. 
    DO t=t1,t2
      Element => Mesh % Elements(t)
      IF( BCMode ) THEN
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE     
        IF ( Element % BoundaryInfo % Constraint /= Tag ) CYCLE
        ActiveNode(Element % NodeIndexes) = .TRUE.
      END IF
    END DO

    ! If all nodes are active just use pointers to the nodes.
    ! Otherwise create list of the nodes. 
    IF( BCMode ) THEN
      NoNodes = COUNT( ActiveNode )
      ALLOCATE( x(NoNodes), y(NoNodes), z(NoNodes) )
      j = 0
      DO i=1,Mesh % NumberOfNodes
        IF(.NOT. ActiveNode(i) ) CYCLE
        j = j + 1
        x(j) = Mesh % Nodes % x(i)
        y(j) = Mesh % Nodes % y(i)
        z(j) = Mesh % Nodes % z(i)
      END DO
    ELSE
      NoNodes = Mesh % NumberOfNodes
      x => Mesh % Nodes % x
      y => Mesh % Nodes % y
      z => Mesh % Nodes % z
    END IF

    ! Call the function to set the sphere parameters for the nodes.
    CALL SphereFitfun(NoNodes,x,y,z,xc,yc,zc,Rad)

    IF( BCMode ) THEN
      DEALLOCATE(ActiveNode,x,y,z)
    END IF

    ! Add the sphere parameters to the list so that they can be used later
    ! directly without having to fit the parameters again.  
    CALL ListAddConstReal( Params,'Sphere Center X',xc )
    CALL ListAddConstReal( Params,'Sphere Center Y',yc )
    CALL ListAddConstReal( Params,'Sphere Center Z',zc )
    CALL ListAddConstReal( Params,'Sphere Radius',Rad )
    
    IF( PRESENT( FitParams ) ) THEN
      FitParams(1) = xc
      FitParams(2) = yc
      FitParams(3) = zc
      FitParams(4) = Rad
    END IF
      
  CONTAINS
    

    ! Sumith YD: Fast Geometric Fit Algorithm for Sphere Using Exact Solution
    !------------------------------------------------------------------------
    SUBROUTINE SphereFitfun(n,x,y,z,xc,yc,zc,R)
      INTEGER :: n
      REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
      REAL(KIND=dp) :: xc,yc,zc,R
      
      REAL(KIND=dp) :: Sx,Sy,Sz,Sxx,Syy,Szz,Sxy,Sxz,Syz,&
          Sxxx,Syyy,Szzz,Syzz,Sxyy,Sxzz,Sxxy,Sxxz,Syyz,&
          A1,a,b,c,d,e,f,g,h,j,k,l,m,delta
      
      Sx = SUM(x); Sy = SUM(y); Sz = SUM(z);
      Sxx = SUM(x*x); Syy = SUM(y*y);
      Szz = SUM(z*z); Sxy = SUM(x*y);
      Sxz = SUM(x*z); Syz = SUM(y*z);
      Sxxx = SUM(x*x*x); Syyy = SUM(y*y*y);
      Szzz = SUM(z*z*z); Sxyy = SUM(x*y*y);
      Sxzz = SUM(x*z*z); Sxxy = SUM(x*x*y);
      Sxxz = SUM(x*x*z); Syyz = SUM(y*y*z);
      Syzz = SUM(y*z*z);

      ! We must do parallel reduction here if the surface is split among
      ! several MPI processes. 
      IF( BCMode .AND. ParEnv % PEs > 1 ) THEN
        Sx = ParallelReduction(Sx); Sy = ParallelReduction(Sy); Sz = ParallelReduction(Sz);
        Sxx = ParallelReduction(Sxx); Syy = ParallelReduction(Syy);
        Szz = ParallelReduction(Szz); Sxy = ParallelReduction(Sxy);
        Sxz = ParallelReduction(Sxz); Syz = ParallelReduction(Syz);
        Sxxx = ParallelReduction(Sxxx); Syyy = ParallelReduction(Syyy);
        Szzz = ParallelReduction(Szzz); Sxyy = ParallelReduction(Sxyy);
        Sxzz = ParallelReduction(Sxzz); Sxxy = ParallelReduction(Sxxy);
        Sxxz = ParallelReduction(Sxxz); Syyz = ParallelReduction(Syyz);
        Syzz = ParallelReduction(Syzz);       
      END IF
           
      A1 = Sxx +Syy +Szz;
      a = 2*Sx*Sx-2*N*Sxx;
      b = 2*Sx*Sy-2*N*Sxy;
      c = 2*Sx*Sz-2*N*Sxz;
      d = -N*(Sxxx +Sxyy +Sxzz)+A1*Sx;
      e = 2*Sx*Sy-2*N*Sxy;
      f = 2*Sy*Sy-2*N*Syy;
      g = 2*Sy*Sz-2*N*Syz;
      h = -N*(Sxxy +Syyy +Syzz)+A1*Sy;
      j = 2*Sx*Sz-2*N*Sxz;
      k = 2*Sy*Sz-2*N*Syz;
      l = 2*Sz*Sz-2*N*Szz;
      m = -N*(Sxxz +Syyz + Szzz)+A1*Sz;
      delta = a*(f*l - g*k)-e*(b*l-c*k) + j*(b*g-c*f);

      xc = (d*(f*l-g*k) -h*(b*l-c*k) +m*(b*g-c*f))/delta;
      yc = (a*(h*l-m*g) -e*(d*l-m*c) +j*(d*g-h*c))/delta;
      zc = (a*(f*m-h*k) -e*(b*m-d*k) +j*(b*h-d*f))/delta;
      R = SQRT(xc*xc+yc*yc+zc*zc+(A1-2*(xc*Sx+yc*Sy+zc*Sz))/N);

    END SUBROUTINE SphereFitfun

  END SUBROUTINE SphereFit


  
  
END MODULE MeshBasics
