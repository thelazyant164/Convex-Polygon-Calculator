# Convex-Polygon-Calculator
My attempt at a program that can receive the Cartesian coordinates of all vertices of a convex polyhedron, then solve for its volume and surface area.

Features:
- The ability to distinguish between convex vertices and concave vertices (those not belonging to the convex hull) in relation to the set of input vertices
- Calculate number of actual edges, faces, vertices of the convex hull created from all input vertices
- Calculate volume and surface area of the convex hull

Limit:
- Can only take manual coordinate input at this stage (can be cumbersome if too many vertices)
- Only limit to 1478 vertices at maximum

Requirements:
- Known all Cartesian coordinates of all vertices

Usage:
- "Convex Polyhedron Solver.rar" is the pre-packaged file meant to run directly on Windows OS. Simply download and extract, then go into "dist" directory and run "Convex Polyhedron Solver.exe"
- "Polyhedron" is the source code, written in Python.
