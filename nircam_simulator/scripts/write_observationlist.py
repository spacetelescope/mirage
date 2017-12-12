# Generate observation list files based on default values and APT output files
from nircam_simulator.scripts import apt_inputs

def write_yaml(xml_file, pointing_file, yaml_file, ps_cat_sw=None, ps_cat_lw=None):
    # Define "default" values (probably should be changed eventually)
    date = '2019-07-04'
    PAV3 = '0.'
    GalaxyCatalog = 'None'
    ExtendedCatalog = 'None'
    ExtendedScale = '1.0'
    ExtendedCenter = '1024,1024'
    MovingTargetList = 'None'
    MovingTargetSersic = 'None'
    MovingTargetExtended = 'None'
    MovingTargetConvolveExtended = 'True'
    MovingTargetToTrack = 'None'
    BackgroundRate_sw = '0.5'
    BackgroundRate_lw = '1.2'

    # Read in filters from XML file
    xml_read = apt_inputs.AptInput()
    xml_table = xml_read.read_xml(xml_file)

    sw_filters = xml_table['ShortFilter']
    lw_filters = xml_table['LongFilter']

    # Read in observation names from pointing file
    f = open(pointing_file)
    rows = f.readlines()
    f.close()
    obs_names = [row[2:-1].split(' (')[0] for row in rows if '* ' == row[:2]]
    num_obs = len(obs_names)

    write = ["# Observation list created by write_observationlist.py script in\n",
             "# nircam_simulator/scripts. Note: all values except filters and\n",
             "# observation names are default.\n\n"]
    for i_obs in range(num_obs):
        write += [\
        "Observation{}:\n".format(i_obs + 1),
        "  Name: '{}'\n".format(obs_names[i_obs]),
        "  Date: {}\n".format(date),
        "  PAV3: {}\n".format(PAV3),
        "  SW:\n",
        "    Filter: {}\n".format(sw_filters[i_obs]),
        "    PointSourceCatalog: {}\n".format(ps_cat_sw),
        "    GalaxyCatalog: {}\n".format(GalaxyCatalog),
        "    ExtendedCatalog: {}\n".format(ExtendedCatalog),
        "    ExtendedScale: {}\n".format(ExtendedScale),
        "    ExtendedCenter: {}\n".format(ExtendedCenter),
        "    MovingTargetList: {}\n".format(MovingTargetList),
        "    MovingTargetSersic: {}\n".format(MovingTargetSersic),
        "    MovingTargetExtended: {}\n".format(MovingTargetExtended),
        "    MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
        "    MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
        "    BackgroundRate: {}\n".format(BackgroundRate_sw),
        "  LW:\n",
        "    Filter: {}\n".format(lw_filters[i_obs]),
        "    PointSourceCatalog: {}\n".format(ps_cat_lw),
        "    GalaxyCatalog: {}\n".format(GalaxyCatalog),
        "    ExtendedCatalog: {}\n".format(ExtendedCatalog),
        "    ExtendedScale: {}\n".format(ExtendedScale),
        "    ExtendedCenter: {}\n".format(ExtendedCenter),
        "    MovingTargetList: {}\n".format(MovingTargetList),
        "    MovingTargetSersic: {}\n".format(MovingTargetSersic),
        "    MovingTargetExtended: {}\n".format(MovingTargetExtended),
        "    MovingTargetConvolveExtended: {}\n".format(MovingTargetConvolveExtended),
        "    MovingTargetToTrack: {}\n".format(MovingTargetToTrack),
        "    BackgroundRate: {}\n\n".format(BackgroundRate_lw)]

    f = open(yaml_file, 'w')
    for line in write:
        f.write(line)
    f.close()
    print('\nSuccessfully wrote {} observations to {}'.format(num_obs, yaml_file))


if __name__ == '__main__':
    xml_file = '../OTECommissioning/OTE01/OTE01-1134.xml'
    pointing_file = '../OTECommissioning/OTE01/OTE01-1134.pointing'
    yaml_file = 'test.yaml'

    write_yaml(xml_file, pointing_file, yaml_file)