import sqlite3
import math as m
import numpy as np
import dynamics as dyn

# Dependencies: tycho.db
#               table name: tycho_data (329499 entries, 325 mb)
#               fields: id, BTmag, VTmag, HIP, RA, DE, name, Bayer,
#               computed_temperature, published_temperature, reduction_term
# Generated Files:

#############################################################
#############################################################
# General Parameters:

# Script options
makeSearchCatalog = True
makeObjectIDCatalog = True
checkNumEntries = True

# BTmag cutoff for generated reference catalogs
M_CUTOFF = 6.

# Decimal places in dtheta entries in objectID catalog
nDecimal = 7

# Maximum angular separation between stars for generated object ID reference catalog
DTHETA_MAX = 15.

#############################################################
#############################################################
# Read original star catalog and generate reference catalogs used in image processing module

if makeSearchCatalog:

    print 'Generating Initial Search Catalog'

    star_catalog = sqlite3.connect('star_catalog/tycho.db')
    with star_catalog:

        s = star_catalog.cursor()

        dec_min = -361
        dec_max = 361

        s.execute("SELECT * FROM tycho_data WHERE BTmag IS NOT NULL "
                  " AND DE BETWEEN (?) AND (?)"
                  " AND BTmag <= (?)",
                  (dec_min, dec_max, M_CUTOFF))
        rows = s.fetchall()

        id = []
        RA = []
        DEC = []
        BTmag = []
        new_rows = []
        for row in rows:
            id.append(row[0])
            RA.append(row[4])
            DEC.append(row[5])
            BTmag.append(row[1])
            new_rows.append((row[0], row[4], row[5], row[1]))


    #######################################################
    # Create new table with subset of original catalog

    objectID_star_catalog = sqlite3.connect('star_catalog/tycho_BTmag_cutoff.db')

    with objectID_star_catalog:

        ob = objectID_star_catalog.cursor()

        table_name = 'tycho_data'

        # delete table if exists
        ob.execute("DROP TABLE IF EXISTS tycho_data")

        # create new table and column for id
        ob.execute('CREATE TABLE {tn} ({nf} {ft} PRIMARY KEY)'
                   .format(tn=table_name, nf='id', ft='INTEGER'))

        # create column for RA
        ob.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
                   .format(tn=table_name, cn='RA', ct='FLOAT'))

        # create column for DEC
        ob.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
                   .format(tn=table_name, cn='DEC', ct='FLOAT'))

        # create column for VTMag
        ob.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
                   .format(tn=table_name, cn='BTMag', ct='FLOAT'))

        for ind_row in range(len(new_rows)):
            ob.execute("INSERT INTO tycho_data(id, RA, DEC, BTmag) VALUES(?,?,?,?)",
                       (new_rows[ind_row]))

    objectID_star_catalog.close()

#######################################
# object ID reference table generation

if makeObjectIDCatalog:

    # pull star catalog entries from reference table
    star_catalog = sqlite3.connect('star_catalog/tycho_BTmag_cutoff.db')
    with star_catalog:

        s = star_catalog.cursor()
        s.execute("SELECT * FROM tycho_data")
        rows = s.fetchall()

        id = []
        ra = []
        dec = []
        ehat = []
        for row in rows:
            id.append(row[0])
            ra.append(row[1])
            dec.append(row[2])

            current_ehat = dyn.radec_to_unitvector((row[1], row[2]))
            ehat.append((current_ehat[0], current_ehat[1], current_ehat[2]))

    # generate angular separation between stars and enter values in stored lists

    id1 = []
    id2 = []
    dtheta = []
    # ra1 = []
    # ra2 = []
    # dec1 =[]
    # dec2 =[]

    n_rows = len(rows)
    n_rows_new = 0

    print 'Calculating Object ID Reference Table'
    # generate row index from the first to second to last
    for ind1 in range(n_rows-1):

        print ind1, ' out of ', n_rows
        # generate row index from current to last index
        for ind2_zero in range(n_rows-ind1-1):

            ind2 = ind2_zero + ind1 + 1

            # calculate angular distance between two row entries
            current_dtheta = m.degrees(m.acos(np.dot(ehat[ind1], ehat[ind2])))

            # store values if angular separation meets criteria
            if current_dtheta <= DTHETA_MAX:

                n_rows_new = n_rows_new + 1
                id1.append(id[ind1])
                id2.append(id[ind2])
                dtheta.append(current_dtheta)
                # ra1.append(ra[ind1])
                # ra2.append(ra[ind2])
                # dec1.append(dec[ind1])
                # dec2.append(dec[ind2])

    id_objectID = range(n_rows_new)


    # create new reference table for object ID purposes
    objectID_catalog = sqlite3.connect('star_catalog/objectID_catalog.db')
    with objectID_catalog:

        print 'Generating Object ID Reference Table'

        oID = objectID_catalog.cursor()

        table_name = 'tycho_objectID'

        # delete table if exists
        oID.execute("DROP TABLE IF EXISTS tycho_objectID")

        # create new table and column for id
        oID.execute('CREATE TABLE {tn} ({nf} {ft} PRIMARY KEY)'
                   .format(tn=table_name, nf='id_objectID', ft='INTEGER'))

        # create column for id1
        oID.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
                   .format(tn=table_name, cn='id1', ct='INTEGER'))

        # create column for id2
        oID.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
                   .format(tn=table_name, cn='id2', ct='INTEGER'))

        # create column for dtheta
        oID.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
                   .format(tn=table_name, cn='dtheta', ct='FLOAT'))

        # create column for RA values
        # oID.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
        #            .format(tn=table_name, cn='RA1', ct='FLOAT'))
        # oID.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
        #            .format(tn=table_name, cn='RA2', ct='FLOAT'))
        #
        # # create column for DEC values
        # oID.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
        #            .format(tn=table_name, cn='DEC1', ct='FLOAT'))
        # oID.execute("ALTER TABLE {tn} ADD COLUMN '{cn}' {ct}"
        #            .format(tn=table_name, cn='DEC2', ct='FLOAT'))

        for ind_row in range(n_rows_new):

            # oID.execute("INSERT INTO tycho_objectID("
            #             "id_objectID, id1, id2, dtheta, RA1, RA2, DEC1, DEC2) VALUES(?,?,?,?,?,?,?,?)",
            #            (id_objectID[ind_row], id1[ind_row], id2[ind_row], dtheta[ind_row],
            #             ra1[ind_row], ra2[ind_row], dec1[ind_row], dec2[ind_row]))
            oID.execute("INSERT INTO tycho_objectID("
                        "id_objectID, id1, id2, dtheta) VALUES(?,?,?,?)",
                       (id_objectID[ind_row], id1[ind_row], id2[ind_row], round(dtheta[ind_row],8)))

    objectID_catalog.close()



if checkNumEntries:

    ref_catalog = sqlite3.connect('star_catalog/objectID_catalog.db')
    with ref_catalog:

        cursor = ref_catalog.cursor()
        cursor.execute("SELECT * FROM tycho_objectID WHERE id1 IS (?) OR id2 IS (?)"
                       "OR id1 IS (?) OR id2 IS (?)", (8291, 8291, 8291, 8291))

        rows = cursor.fetchall()

        print 'Number of Entries: ', len(rows)
        for row in rows:
            print row

    ref_catalog = sqlite3.connect('star_catalog/tycho.db')
    with ref_catalog:

        cursor = ref_catalog.cursor()
        cursor.execute("SELECT * FROM tycho_data WHERE BTmag IS NOT NULL "
                  " AND BTmag <= (?)", (M_CUTOFF,))
        rows = cursor.fetchall()

        nStars = len(rows)
        averagePerSqDeg = nStars / 41253.

        print 'Number of Catalog Entries for BTmag <= ', M_CUTOFF, ': ', nStars
        print 'Average # of Stars Per 100 sq deg: ', averagePerSqDeg*100


