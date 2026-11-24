""" Tests related to FRB ingestion (e.g. the mw_ebv IRSA dust query
made in YSE_App.models.frbtransient_models.execute_after_save) """

from astropy.coordinates import SkyCoord
from astroquery import irsa_dust


def get_mw_ebv(ra: float = 150.0, dec: float = 30.0) -> float:
    """ Recover the Galactic E(B-V) at (ra, dec) from the IRSA dust service,
    using the same call as execute_after_save()

    Args:
        ra (float): RA in deg
        dec (float): Dec in deg

    Returns:
        float: mw_ebv (Schlafly & Finkbeiner mean)
    """
    c = SkyCoord(ra=ra, dec=dec, unit='deg')
    tbl = irsa_dust.IrsaDust.get_query_table(c)
    mw_ebv = tbl['ext SandF mean'][0]
    return float(mw_ebv)


def test_mw_ebv():
    mw_ebv = get_mw_ebv()
    assert 0. <= mw_ebv < 1.


if __name__ == '__main__':
    # Diagnostics -- useful when run inside the web container:
    #   docker exec ysepz_web_container python3 testing/test_ingest.py
    import sys
    import certifi
    import astroquery
    print(f"python {sys.version.split()[0]}, astroquery {astroquery.__version__}, "
          f"certifi {certifi.__version__}")
    print(f"mw_ebv = {get_mw_ebv()}")
