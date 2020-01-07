# Geospatial Function Files
# Author: Nicole Hemming-Schroeder
# Email: hemmingn at uci.edu

def albers2latlon(x, y, silent = 1):
    ## Algorithm from John Parr Synder, USGS, 1987, for the ellipsoid
    # Publication title: Map Projections -- a working manual: issue 1395
    # Free from Google Play: https://play.google.com/store/books/details?id=numpydOAAAAMAAJ&rdid=book-numpydOAAAAMAAJ&rdot=1

    # Inputs: (x, y, silent)
    # x: Albers x-coordinate
    # y: Albers y-coordinate
    # silent: default (1), set to 0 if you want the code to tell you whether the longitude estimates converged or not 
    
    # Outputs: latitude, longitude
    
    # Checking the code against Synder example (see table values on pg. 103)
    # x=1885472.7
    # y=1535925.0

    import numpy
    # The following parameters are all consistent with the Albers projection described in the Wilson file as well
    # as the USGS Synder paper (1987)
    phi_0 = 23 # degrees
    lambda_0 = -96 # degrees
    phi_1 = 29.5 # degrees
    phi_2 = 45.5 # degrees
    a = 6378206.4 # radius of the Earth's semi-major axis (Equator) in meters 
    ee = 0.0822719 # ? (see pg. 292 of Synder (1987) - give it a double ee to avoid confusion with the standard e)

    m_1 = numpy.cos(phi_1*numpy.pi/180)/numpy.sqrt(1-ee**2*numpy.sin(phi_1*numpy.pi/180)**2)
    m_2 = numpy.cos(phi_2*numpy.pi/180)/numpy.sqrt(1-ee**2*numpy.sin(phi_2*numpy.pi/180)**2)

    # note numpy.log refers to the natural logarithm
    q_0 = (1-ee**2)*numpy.abs(numpy.sin(phi_0*numpy.pi/180)/(1-(ee**2)*numpy.sin(phi_0*numpy.pi/180)**2) - (1/(2*ee))*numpy.log((1-               ee*numpy.sin(phi_0*numpy.pi/180))/(1+ee*numpy.sin(phi_0*numpy.pi/180))))
    q_1 = (1-ee**2)*numpy.abs(numpy.sin(phi_1*numpy.pi/180)/(1-(ee**2)*numpy.sin(phi_1*numpy.pi/180)**2) - (1/(2*ee))*numpy.log((1-             ee*numpy.sin(phi_1*numpy.pi/180))/(1+ee*numpy.sin(phi_1*numpy.pi/180))))
    q_2 = (1-ee**2)*numpy.abs(numpy.sin(phi_2*numpy.pi/180)/(1-(ee**2)*numpy.sin(phi_2*numpy.pi/180)**2) - (1/(2*ee))*numpy.log((1-ee*numpy.sin(phi_2*numpy.pi/180))/(1+ee*numpy.sin(phi_2*numpy.pi/180))))

    # Compute C and n from the previous
    n  = (m_1**2 - m_2**2)/(q_2 - q_1)
    C = m_1**2 + n*q_1
    rho_0 = a*numpy.sqrt(C - n * q_0)/n
    rho = numpy.sqrt(x**2 + (rho_0 - y)**2)
    theta = numpy.arctan(x/(rho_0-y))*180/numpy.pi 
    q = (C - rho**2*n**2/a**2)/n

    # Compute latitude (phi) and longitude (lambda)
    # Need to use iterative equation (3-16) to find phi

    # Initial guess:
    phi = numpy.arcsin(q/2)*180/numpy.pi
    
    # Iterate
    nits = 100
    for i in numpy.arange(0,nits):
        #print(phi)
        phi_next = phi + ((1 - ee**2*numpy.sin(phi*numpy.pi/180)**2)**2/(2*numpy.cos(phi*numpy.pi/180)))*(q/(1-ee**2) - numpy.sin(phi*numpy.pi/180)/(1-ee**2*numpy.sin(phi*numpy.pi/180)**2) + 1/(2*ee)*numpy.log((1-ee*numpy.sin(phi*numpy.pi/180))/(1+ee*numpy.sin(phi*numpy.pi/180))))*180/numpy.pi

        error = phi - phi_next
        # Impose |error| less than 1e-7
        if abs(error)<1e-7:
            if silent == 0:
                print('converged in ' + str(i) + ' iterations')
            lambda_lon = (lambda_0 + theta/n)
            return(phi, lambda_lon)
            break
        phi = phi_next*1.0
        
    if numpy.abs(error)>=1e-7:
        if silent == 0:
            print('did not converge')
        lambda_lon = numpy.nan
        return(phi, lambda_lon)
        

def wgs84_radius(phi):
    # Finds the distance from the center of the WGS 84 ellipsoid to the surface of the
    # Earth for a given latitude phi, given in degrees
    # Returns the distance in meters
    import numpy
    a = 6378137 # major axis at Equator in meters
    b = 6356752.3142 # minor axis at Poles in meters
    d = numpy.sqrt((a**2)*(b**2)/(a**2 + (b**2 - a**2)*numpy.cos(phi*numpy.pi/180)**2))
    return(d)
