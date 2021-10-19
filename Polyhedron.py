#Import libraries and functions
from itertools import combinations
from math import sqrt
from copy import deepcopy
from numpy.linalg import solve, LinAlgError
from fractions import Fraction



#Define algebra toolbox

#Define method to solve system of two linear equations
#left = [row1, row2]
#right = [numRow1, numRow2]
#row1 = [codependent1, codependent2]
#X = [x, y]

#Using numpy.linalg.solve()
def sysLinearTwo(left, right):
    A = left
    B = right
    X = solve(A, B)
    return X



#Define geometry toolbox

#Define function to create unit vector between two given points
#unitVector = [x2 - x1, y2 - y1, z2 - z1]
def unitVector(point1, point2):
    unitVector = []
    for i in range(3):
        unitVector.append(point2[i] - point1[i])
    return unitVector

#Define function to create normal vector for two given vectors
#normalVector = [b1c2 - c1b2, c1a2 - a1c2, a1b2 - a2b1]
def normalVector(vector1, vector2):
    a = vector1[1]*vector2[2] - vector1[2]*vector2[1]
    b = vector1[2]*vector2[0] - vector1[0]*vector2[2]
    c = vector1[0]*vector2[1] - vector1[1]*vector2[0]
    normalVector = [a, b, c]
    return normalVector

#Define function to determine linear equation containing two given points
#line = [unitVector, point1]
def linear(point1, point2):
    line = []
    unit = unitVector(point1, point2)
    line.extend((unit, point1))
    return line

#Define function to check if a point belongs in a given line
def pointInLine(pointToCheck, line):
    verifier = []
    coefficientTracker = []
    for i in range(len(line[1])):
        if (line[0][i] == 0):
            if (pointToCheck[i] == line[1][i]):
                verifier.append(True)
            else:
                verifier.append(False)
        else:
            coefficientTracker.append(Fraction((pointToCheck[i] - line[1][i])/line[0][i]).limit_denominator(1000000))
            verifier.append(True)            
    if (len(coefficientTracker) > 1):
        for coefficient in coefficientTracker:
            if (coefficientTracker[0] != coefficient): return False
    return (verifier[0] and verifier[1] and verifier[2])

#Define function to determine plane equation containing three given points
#plane = [normalVector, codependent]
#normalVector = [x, y, z]
def plane(point1, point2, point3):
    
    #Making certain that the three given points don't align
    if (pointInLine(point1, linear(point2, point3))):
        exit

    plane = []
    unitVector1 = unitVector(point1, point2)
    unitVector2 = unitVector(point1, point3)
    normal = normalVector(unitVector1, unitVector2)
    codependent = - (point1[0] * normal[0] + point1[1] * normal[1] + point1[2] * normal[2])
    plane.extend((normal, codependent))
    return plane

#Define function to check if two given planes are parallel
def isParallel(plane1, plane2):
    verifier = []
    coefficientTracker = []
    for i in range(len(plane1[0])):
        if (plane2[0][i] == 0):
            if (plane1[0][i] == 0):
                verifier.append(True)
            else:
                verifier.append(False)
        else:
            coefficientTracker.append(Fraction(plane1[0][i]/plane2[0][i]).limit_denominator(1000000))
            verifier.append(True)
    if (len(coefficientTracker) > 1):
        for coefficient in coefficientTracker:
            if (coefficientTracker[0] != coefficient): return False
    return (verifier[0] and verifier[1] and verifier[2])

#Define function to check if two given planes are identical
def isIdentical(plane1, plane2, inputCoords):
    if (not isParallel(plane1, plane2)):
        return False
    for coord in inputCoords:
        if (pointInPlane(coord, plane1)):
            point = coord
            break
    return (pointInPlane(point, plane2))

#Define function to determine linear equation of two planes' intersection
def intersection(plane1, plane2):

    #Determine normal vector for intersection
    normal1 = plane1[0]
    normal2 = plane2[0]
    normal = normalVector(normal1, normal2)

    #Determine a point which intersection travels through
    while (True):
    #Check for special cases
    #Case 1: intersection travels through a point with coordinates (x, 0, z)
        a1 = float(plane1[0][0])
        a2 = float(plane2[0][0])
        b1 = float(plane1[0][2])
        b2 = float(plane2[0][2])
        c1 = float(- plane1[1])
        c2 = float(- plane2[1])
        left = [[a1, b1], [a2, b2]]
        right = [c1, c2]
        try:
            point = sysLinearTwo(left, right).tolist()
            point.insert(1, 0)
            break
        except LinAlgError:
            pass
    #Case 2: intersection travels through a point with coordinates (0, y, z)
        a1 = float(plane1[0][1])
        a2 = float(plane2[0][1])
        b1 = float(plane1[0][2])
        b2 = float(plane2[0][2])
        c1 = float(- plane1[1])
        c2 = float(- plane2[1])
        left = [[a1, b1], [a2, b2]]
        right = [c1, c2]
        try:
            point = sysLinearTwo(left, right).tolist()
            point.insert(0, 0)
            break
        except LinAlgError:
            pass
    #Case 3: intersection travels through a point with coordinates (x, y, 0)
        a1 = float(plane1[0][0])
        a2 = float(plane2[0][0])
        b1 = float(plane1[0][1])
        b2 = float(plane2[0][1])
        c1 = float(- plane1[1])
        c2 = float(- plane2[1])
        left = [[a1, b1], [a2, b2]]
        right = [c1, c2]
        try:
            point = sysLinearTwo(left, right).tolist()
            point.insert(2, 0)
            break
        except LinAlgError:
            pass

    #A line cannot be parallel to all three planes x = 0, y = 0, z = 0 at the same time, so at least one of the three special cases yielded viable products
    intersection = []

    #Revert numpy-array's float-type return back into fraction-type for higher precision
    point[0] = Fraction(point[0]).limit_denominator(1000000)
    point[1] = Fraction(point[1]).limit_denominator(1000000)
    point[2] = Fraction(point[2]).limit_denominator(1000000)
    intersection.extend((normal, point))
    return intersection

#Define function to find a point's relative position to a plane in numerical value
def relativePosition(point, plane):
    normal = plane[0]
    codependent = plane[1]
    relativePosition = normal[0] * point[0] + normal[1] * point[1] + normal[2] * point[2] + codependent
    return relativePosition

#Define function to check if a given point belongs in a given plane
def pointInPlane(point, plane):
    return (relativePosition(point, plane) == 0)

#Define function to find distance between two given points
def distance(point1, point2):
    unit = unitVector(point1, point2)
    distance = sqrt(unit[0]**2 + unit[1]**2 + unit[2]**2)
    return distance

#Define function to calculate distance between given point and plane
def calculateDistance(point, plane):
    normal = plane[0]
    x = point[0]
    y = point[1]
    z = point[2]
    a = normal[0]
    b = normal[1]
    c = normal[2]
    d = plane[1]
    distance = Fraction(abs(a*x + b*y + c*z + d)/sqrt(a**2 + b**2 + c**2)).limit_denominator(1000000)
    return distance



#Define program-specific functions

#Define function to verify and receive input of all vertices' coordinates and append them to inputCoords
#inputCoords = [point1, point2...]
#point1 = [x, y, z]
def inputVertexCoords():
    inputCoords = []
    valid = False
    while (valid == False):
        count = float(input("How many vertices does this polyhedron have? "))

        #Check to see if vertices count is valid
        if (count > 1477 or count < 4 or count - int(count) != 0):
            valid = False
            print('Invalid vertex count. Only an integer value between 4 and 1477 may be accepted.\nPlease try again.\n')
        else:
            valid = True

    first = True
    i = 0
    while i < int(count):
        newCoord = []

        #Receive input and append to list as a sublist
        x, y, z = input(f"Input coordinates x, y and z of point #{i + 1} respectively, separated by a comma: ").split(',')
        newCoord.extend((Fraction(x), Fraction(y), Fraction(z)))
        i += 1

        #Pass checking for the first time
        if (first):
            inputCoords.append(newCoord)
            first = False
            continue

        #Check to see if any duplicated vertices were input
        try:
            inputCoords.index(newCoord)
            print("Duplicated vertices detected. Please input another vertex.")
            i -= 1
            continue
        except ValueError:
            inputCoords.append(newCoord)
            continue

    return inputCoords

#Define function to list all possible plane equations from inputCoords
def generateAllPlanes(inputCoords):
    allPlanes = []
    trios = combinations(inputCoords, 3)
    for trio in trios:
        allPlanes.append(plane(trio[0], trio[1], trio[2]))
    return [plane for plane in allPlanes if (plane is not None and plane != [[0, 0, 0], 0])]

#Define function to iterate through all inputCoords and add all planes containing faces to list:
def findFace(inputCoords):

    #Define function to check if a plane contains one of the polyhedron's faces - by being on one side compared to all other vertices
    #points = inputCoords
    def containsFace(points, plane):

        #Create mirror list containing all points not belonging in plane
        outsidePoints = []
        for point in points:
            if (pointInPlane(point, plane)):
                continue
            outsidePoints.append(point)

        #Check the multiplication of each point's relative position to plane with the first point's relative position to plane
        for point in outsidePoints:
            if (relativePosition(point, plane) * relativePosition(outsidePoints[0], plane) < 0): return False
        return True

    allPlanes = generateAllPlanes(inputCoords)
    facePlanes = []
    for plane in allPlanes:
        if (containsFace(inputCoords, plane)):
            facePlanes.append(plane)

    #Making certain each face in facePlanes is unique
    duos = combinations(facePlanes, 2)
    for duo in duos:
        if (isIdentical(duo[0], duo[1], inputCoords)):
            try:
                facePlanes.index(duo[0])
                facePlanes.remove(duo[1])
            except ValueError:
                pass
    return facePlanes

#Define function to verify actual vertices (belonging to three faces), invalid (belonging to two or one faces) or eliminated (belonging to none)
def verifyActualVertices(inputCoords, facePlanes):
    actualVertices = []
    invalidVertices = []
    eliminatedVertices = []
    for vertex in inputCoords:
        counter = 0
        for face in facePlanes:
            if (pointInPlane(vertex, face)): counter += 1
        if (counter >= 3):
            actualVertices.append(vertex)
        elif (counter >= 1):
            invalidVertices.append(vertex)
        else:
            eliminatedVertices.append(vertex)
    return actualVertices, invalidVertices, eliminatedVertices

#Define function to list all possible plane intersections from facePlanes
def generateAllIntersections(facePlanes):
    allIntersections = []
    duos = combinations(facePlanes, 2)
    for duo in duos:
        if (not isParallel(duo[0], duo[1])):
            allIntersections.append(intersection(duo[0], duo[1]))
    return allIntersections

#Define function to identify and add all edges to list
def findEdges(actualVertices, facePlanes):

    #Define function to check if an intersection goes through at least two vertices from actualVertices
    def verifyEdge(intersection, actualVertices):
        counter = 0
        for point in actualVertices:
            if (pointInLine(point, intersection)): counter += 1
        return (counter >= 2)

    foundEdges = []
    allIntersections = generateAllIntersections(facePlanes)
    for intersection in allIntersections:
        if (verifyEdge(intersection, actualVertices)): foundEdges.append(intersection)
    return foundEdges

#Define optional function to display all invalid vertices
def showInvalidVertices(invalidVertices, inputCoords):
    if (len(invalidVertices) != 0):
        choice = input(f'\n{len(invalidVertices)} invalid vertices detected. Do you want to view them? (y) for yes or (n) for no: ')
        if (choice == 'n'):
            exit
        for invalidVertex in invalidVertices:
            displayInvalidVertex = []
            for coord in invalidVertex:
                displayInvalidVertex.append(float(coord))
            print(f'Point #{inputCoords.index(invalidVertex) + 1}, with coordinates {displayInvalidVertex}, was invalid. This vertex belonged to a pre-existing edge or face, and thus was disqualified from being computed.\n')
        exit
    exit

def showEliminatedVertices(eliminatedVertices, inputCoords):
    if (len(eliminatedVertices) != 0):
        choice = input(f'\n{len(eliminatedVertices)} eliminated vertices detected. Do you want to view them? (y) for yes or (n) for no: ')
        if (choice == 'n'):
            exit
        for eliminatedVertex in eliminatedVertices:
            displayEliminatedVertex = []
            for coord in eliminatedVertex:
                displayEliminatedVertex.append(float(coord))
            print(f'Point #{inputCoords.index(eliminatedVertex) + 1}, with coordinates {displayEliminatedVertex}, was eliminated. This vertex belonged inside the convex hull of the polyhedron created by the given coordinates, and thus was disqualified from being computed.\n')
        exit
    exit

#Define function to verify if actual vertices connect to create a convex polyhedron
def verifyConvexPolyhedron(actualVertices, facePlanes):

    #Define function to verify that a trihedral vertex exists
    def verifyTrihedral(actualVertices, facePlanes):
        for vertex in actualVertices:
            adjacentFaceCounter = 0
            for face in facePlanes:
                if (pointInPlane(vertex, face)): adjacentFaceCounter += 1
            if (adjacentFaceCounter == 3): return True
        return False
    
    #Define function to verify that a triangular face exists
    def verifyTriangular(actualVertices, facePlanes):
        for face in facePlanes:
            containedVertexCounter = 0
            for vertex in actualVertices:
                if (pointInPlane(vertex, face)): containedVertexCounter += 1
            if (containedVertexCounter == 3): return True
        return False
    
    if (verifyTriangular(actualVertices, facePlanes) or verifyTrihedral(actualVertices, facePlanes)):
        print("\nInput coordinates allow for possible convex polyhedron creation.\n")
        return True
    else:
        print("\nInput coordinates do not allow for possible convex polyhedron creation.")
        print('Calculations failed: the input coordinates cannot form a plausible convex polyhedron.')
        return False

#Define function to calculate all faces' surface area and add to list
def calculateSurfaceArea(actualVertices, facePlanes):
    
    #Define function to split all vertices into corresponding group with each face's plane
    def groupFaces(actualVertices, facePlanes):
        groupedVertices = []
        for face in facePlanes:
            group = []
            vertexGroup = []
            for vertex in actualVertices:
                if (pointInPlane(vertex, face)): vertexGroup.append(vertex)
            group.extend((vertexGroup, face))
            groupedVertices.append(group)
        return groupedVertices
    
    #Define function to calculate area of triangle from three given points
    def calculateTriangleArea(point1, point2, point3):
        unit1 = unitVector(point1, point2)
        unit2 = unitVector(point1, point3)
        normal = normalVector(unit1, unit2)
        triangleArea = Fraction(sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)/2).limit_denominator(1000000)
        return triangleArea

    #Define function to calculate surface area of polygon on one plane from given vertices
    def calculateFlatArea(vertexGroup):
        vertexGroupCopy = deepcopy(vertexGroup)
        flatArea = 0
        arrangedVertexGroup = [vertexGroupCopy.pop(0)]
        while (len(arrangedVertexGroup) < len(vertexGroup)):
            vertexDistance = []
            for vertex in vertexGroupCopy:
                vertexDistance.append(distance(arrangedVertexGroup[-1], vertex))
            for i in range(len(vertexDistance)):
                if (vertexDistance[i] == min(vertexDistance)):
                    arrangedVertexGroup.append(vertexGroupCopy.pop(i))
                    break
        for i in range(len(arrangedVertexGroup) - 2):
            flatArea += calculateTriangleArea(arrangedVertexGroup[0], arrangedVertexGroup[i + 1], arrangedVertexGroup[i + 2])
        return flatArea

    #Return list of separate faces' surface area for volume calculation
    facesArea = []
    groupedFaces = groupFaces(actualVertices, facePlanes)
    for groupedFace in groupedFaces:
        areaFace = []
        vertices = groupedFace[0]
        facePlane = groupedFace[1]
        areaFace.extend((calculateFlatArea(vertices), facePlane))
        facesArea.append(areaFace)

    #Return actual surface area, totaled
    totalArea = 0
    for faceArea in facesArea:
        totalArea += faceArea[0]
    return facesArea, totalArea

#Define function to calculate total volume by splitting polyhedron into multiple pyramidal polyhedrons from a given vertex
def calculateTotalVolume(facesArea):

    #Define function to calculate a given polygonal pyramid's volume
    def calculateVolume(faceArea, point):
        volume = Fraction(1/3 * faceArea[0] * calculateDistance(point, faceArea[1])).limit_denominator(1000000)
        return volume

    #Define function to pick a random point inside polyhedron
    def generatePointInPolyhedron(actualVertices):
        verticesDistance = []
        duos = combinations(actualVertices, 2)
        for duo in duos:
            verticesDistance.append(distance(duo[0], duo[1]))
        for i in range(len(verticesDistance)):
            if (verticesDistance[i] == max(verticesDistance)):
                point1 = duo[0]
                point2 = duo[1]
        x = Fraction((point1[0] + point2[0])/2).limit_denominator(1000000)
        y = Fraction((point1[1] + point2[1])/2).limit_denominator(1000000)
        z = Fraction((point1[2] + point2[2])/2).limit_denominator(1000000)
        midpoint = [x, y, z]
        return midpoint

    point = generatePointInPolyhedron(actualVertices)
    totalVolume = 0
    for face in facesArea:
        totalVolume += calculateVolume(face, point)
    return totalVolume



#Program flow
inputCoords = inputVertexCoords()
facePlanes = findFace(inputCoords)
actualVertices, invalidVertices, eliminatedVertices = verifyActualVertices(inputCoords, facePlanes)
foundEdges = findEdges(actualVertices, facePlanes)
checkpoint = verifyConvexPolyhedron(actualVertices, facePlanes)
if (checkpoint):
    facesArea, totalArea = calculateSurfaceArea(actualVertices, facePlanes)
    totalVolume = calculateTotalVolume(facesArea)

    #Output
    print(f'Polyhedron has {len(facePlanes)} faces.')
    print(f'Polyhedron has {len(foundEdges)} edges.')
    print(f'Polyhedron has {len(actualVertices)} actual vertices.')
    showInvalidVertices(invalidVertices, inputCoords)
    showEliminatedVertices(eliminatedVertices, inputCoords)
    answer = input("Do you wish to view the results in decimal form or fraction form?\n(f) for fraction or (d) for decimal: ")
    if (answer == "d"):
        accuracy = int(input("How high do you want the calculated precision to be?\n(Enter the number of decimal value to calculate accurately towards, maximum 10000000) "))
        print(f'Polyhedron has a total surface area of {round(float(totalArea), accuracy)}.')
        print(f'Polyhedron has a total volume of {round(float(totalVolume), accuracy)}.')
    else:
        print(f'Polyhedron has a total surface area of {(totalArea)}.')
        print(f'Polyhedron has a total volume of {(totalVolume)}.')