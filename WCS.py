import numpy as np



def WCS(header, shape):
    '''
    read the keywords in header which are used to calculate the wcs and do the 3 steps to find the euqarotial coordinate for the pixels in the cube
    :param header: cube's header
    :param shape: cube's shape
    :return: equarotial wcs
    '''

    # read header
    referpoint_physical = np.array([header['CRPIX1'], header['CRPIX2']])
    referpoint_equarotial = np.array([header['CRVAL1'], header['CRVAL2']])
    phi_p = header['LONPOLE']*np.pi/180.
    transfer_matrix = np.array([[header['PC1_1'], header['PC1_2']], [
                               header['PC2_1'], header['PC2_2']]])

    # create a n*m*2 array to hold coordinate for every pixel
    wcs_map = np.zeros((shape[2], shape[1], 2))
    # put the physical coordinate of the pixel in corresponding cell which will be used to generate intermediate wcs
    for i in range(shape[2]):
        for j in range(shape[1]):
            wcs_map[i, j] = [i+1, j+1]

    # calculate intermediate wcs
    intermediate_wcs = Physical2Intermediate(
        wcs_map, transfer_matrix, referpoint_physical)
    # calculate sphere wcs with intermediate wcs
    sphere_wcs = DZenithalProjection(intermediate_wcs)
    # calculate equarotial wcs with sphere wcs
    equarotial_wcs = CoordinateRotation(
        sphere_wcs, referpoint_equarotial, phi_p)

    return equarotial_wcs


def Physical2Intermediate(wcs_map, transfer_matrix, referpoint_physical):
    '''
    convert physical coordinate to intermediate coordinates by transformation matrix in the header
    :param wcs_map: 3D array contain physical coordinate
    :param transfer_matrix: this matrix used to convert physical coordinate to intermediate coordinate
    :param referpoint_physical: physical coordinate of reference point
    :return: intermediate wcs
    '''
    deta_wcs_map = wcs_map-referpoint_physical
    intermediate_wcs = np.dot(deta_wcs_map, transfer_matrix.T)
    return intermediate_wcs


def DZenithalProjection(intermediate_wcs):
    '''
    convert intermediate coordinate to sphere coordinate
    natue of this step is porjecting 2D plane coordinate system to 2D sphere coordinate system
    :param intermediate_wcs:
    :return: sphere coordinate system
    '''
    intermediate_wcs_shape = np.shape(intermediate_wcs)
    sphere_wcs = np.zeros(intermediate_wcs_shape)

    # do the convertion for every pixels
    # this is zenithal projection which is a very simple projection(rules of projection can be found in wiki)
    for i in range(intermediate_wcs_shape[0]):
        for j in range(intermediate_wcs_shape[1]):
            r_theta = np.sqrt(np.sum(intermediate_wcs[i, j]**2))
            sphere_wcs[i, j] = [np.arctan2(-intermediate_wcs[i, j, 0], intermediate_wcs[i, j, 1]),
                                np.arctan(180. / (r_theta*np.pi))]
            # sphere_wcs[i, j] = [np.arctan2(-intermediate_wcs[i,j,0],intermediate_wcs[i,j,1]), np.arctan(180. / (r_theta * np.pi))]
    # sphere_wcs=sphere_wcs*180./np.pi
    return sphere_wcs


def CoordinateRotation(sphere_wcs, referpoint_equarotial, phi_p):
    '''
    this step rotates the axe of the sphere to the posstion of axes of equarotial coordinate system
    :param sphere_wcs:
    :param referpoint_equarotial: equarotial coordinate of reference point
    :param phi_p:
    :return: equarotial coordinate
    '''

    referpoint_equarotial = referpoint_equarotial*np.pi/180.
    sphere_wcs_shape = np.shape(sphere_wcs)
    equarotial_wcs = np.zeros(sphere_wcs_shape)

    for i in range(sphere_wcs_shape[0]):
        for j in range(sphere_wcs_shape[1]):
            tan = (np.cos(sphere_wcs[i, j, 1])*np.sin(sphere_wcs[i, j, 0]))/(np.sin(sphere_wcs[i, j, 1])*np.cos(
                referpoint_equarotial[1])+np.cos(sphere_wcs[i, j, 1])*np.cos(sphere_wcs[i, j, 0])*np.sin(referpoint_equarotial[1]))
            sin = np.sin(sphere_wcs[i, j, 1])*np.sin(referpoint_equarotial[1])+np.cos(
                sphere_wcs[i, j, 1])*np.cos(referpoint_equarotial[1])*np.cos(sphere_wcs[i, j, 0]-phi_p)
            equarotial_wcs[i, j] = [
                referpoint_equarotial[0]+np.arctan(tan), np.arcsin(sin)]
    equarotial_wcs[13, 47] = referpoint_equarotial

    return equarotial_wcs*180./np.pi