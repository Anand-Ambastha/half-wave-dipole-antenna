' ============================================================
' Half-wave dipole at 2.4 GHz: CST Microwave Studio macro
' Run: CST MWS > Home > Macros > (or File > Macro > Run Macro)
' All dimensions in mm, frequency in GHz.
' Parameters are stored in the Parameter List so the Parameter
' Sweep / Optimizer can vary them (L is the main tuning knob).
' ============================================================
Sub Main ()

'--- units and frequency range ---------------------------------
With Units
    .Geometry "mm"
    .Frequency "GHz"
    .Time "ns"
End With
Solver.FrequencyRange "1.5", "3.5"

'--- parameters ------------------------------------------------
StoreDoubleParameter "f0", 2.4
StoreDoubleParameter "L",  57.1      ' tip-to-tip length (start value from Python MoM resonance)
StoreDoubleParameter "a",  1.0       ' wire radius
StoreDoubleParameter "g",  2.0       ' feed gap
RebuildOnParametricChange False, False

'--- background and boundaries ---------------------------------
With Background
    .Type "Normal"
    .Epsilon "1.0"
    .Mu "1.0"
End With
With Boundary
    .Xmin "expanded open" : .Xmax "expanded open"
    .Ymin "expanded open" : .Ymax "expanded open"
    .Zmin "expanded open" : .Zmax "expanded open"
    .Xsymmetry "none" : .Ysymmetry "none" : .Zsymmetry "none"
End With

'--- material: PEC (switch to Copper (annealed) for loss study) -
' (PEC is built in; nothing to define)

'--- upper arm -------------------------------------------------
With Cylinder
    .Reset
    .Name "arm_top"
    .Component "dipole"
    .Material "PEC"
    .OuterRadius "a"
    .InnerRadius "0"
    .Axis "z"
    .Zrange "g/2", "L/2"
    .Xcenter "0"
    .Ycenter "0"
    .Segments "0"
    .Create
End With

'--- lower arm -------------------------------------------------
With Cylinder
    .Reset
    .Name "arm_bottom"
    .Component "dipole"
    .Material "PEC"
    .OuterRadius "a"
    .InnerRadius "0"
    .Axis "z"
    .Zrange "-L/2", "-g/2"
    .Xcenter "0"
    .Ycenter "0"
    .Segments "0"
    .Create
End With

'--- feed: discrete port across the gap, 50 ohm -----------------
With DiscretePort
    .Reset
    .PortNumber "1"
    .Type "SParameter"
    .Impedance "50.0"
    .Voltage "1.0"
    .Current "1.0"
    .Monitor "True"
    .Radius "0.0"
    .SetP1 "False", "0", "0", "-g/2"
    .SetP2 "False", "0", "0", "g/2"
    .InvertDirection "False"
    .LocalCoordinates "False"
    .Wire ""
    .Position "end1"
    .Create
End With

'--- far-field monitors ----------------------------------------
Dim fr As Variant
For Each fr In Array(2.2, 2.4, 2.6)
    With Monitor
        .Reset
        .Name "farfield (f=" & CStr(fr) & ")"
        .Domain "Frequency"
        .FieldType "Farfield"
        .Frequency CStr(fr)
        .UseSubvolume "False"
        .Create
    End With
Next fr

'--- mesh ------------------------------------------------------
With Mesh
    .MeshType "PBA"
    .SetCreator "High Frequency"
End With
With MeshSettings
    .SetMeshType "Hex"
    .Set "Version", 1%
    .Set "StepsPerWaveNear", "20"
    .Set "StepsPerWaveFar",  "20"
    .Set "StepsPerBoxNear",  "20"
    .Set "StepsPerBoxFar",   "20"
End With
With MeshSettings
    .SetMeshType "Hex"
    .Set "RatioLimitGeometry", "20"
End With

'--- solver (time domain, -40 dB accuracy) ----------------------
ChangeSolverType "HF Time Domain"
With Solver
    .Method "Hexahedral"
    .StimulationPort "All"
    .StimulationMode "All"
    .SteadyStateLimit "-40"
    .MeshAdaption False
    .CalculateModesOnly False
End With

Rebuild
Solver.Start

End Sub
