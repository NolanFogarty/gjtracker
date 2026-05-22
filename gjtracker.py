# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "flask>=3.0",
# ]
# ///
"""GJTRACKER - Moon/Sun tracker.

Python port of GJTRACKER.EXE (Visual Basic 6) by Lance Collister, W7GJ.
Provides a Flask web GUI for computing moon/sun position tracking schedules
for EME (Earth-Moon-Earth) ham radio operation.

Run with:  uv run gjtracker.py
"""

import math
from datetime import date, timedelta
from flask import Flask, request, jsonify, render_template_string


# =============================================================================
# Constants
# =============================================================================
PI = math.pi
TUPI = 2 * PI
RAD = TUPI / 360       # degrees -> radians
DEG = 360 / TUPI       # radians -> degrees
EARTH_RADIUS_KM = 6378.16


# =============================================================================
# Math helpers (preserved names from VB)
# =============================================================================
def fna(x_rad):
    """Radians -> degrees, rounded to nearest tenth (VB FNA)."""
    return int(x_rad * DEG * 10 + 0.5) / 10


def fnc(x):
    """Fractional turns -> radians in [0, 2pi) (VB FNC)."""
    return (x - int(x)) * TUPI


def fnjulian(ay, am, ad):
    """Julian date minus 2397547.5 for 0000 UTC (VB FNJULIAN, years 1900-2099)."""
    return 367 * ay - int(7 * (ay + int((am + 9) / 12)) / 4) + int(275 * am / 9) + ad - 676534


def fnpol(az, el, lat):
    """Cross polarization angle wrt earth axis in degrees (VB FNPOL)."""
    return int(math.atan(
        (math.sin(lat) * math.cos(el) - math.cos(lat) * math.cos(az) * math.sin(el))
        / (math.cos(lat) * math.sin(az))
    ) * DEG)


# =============================================================================
# Maidenhead locator <-> lat/lon
# =============================================================================
def locator_to_latlon(locator):
    """Convert Maidenhead locator (e.g. 'CN87xq') to (lat_deg, lon_deg).

    Returns lat (positive N, negative S) and lon (positive E, negative W).
    """
    loc = locator.strip().upper()
    if len(loc) == 4:
        loc += "MM"
        x = 0
    elif len(loc) == 2:
        loc += "55MM"
        x = 0
    elif len(loc) >= 6:
        loc = loc[:6]
        x = 0.5
    else:
        raise ValueError("Locator must be 2, 4, or 6 characters")

    lon_pos_w = 180 - (ord(loc[0]) - 65) * 20 - int(loc[2]) * 2 - ((ord(loc[4]) - 65) + x) / 12
    lat_pos_n = -90 + (ord(loc[1]) - 65) * 10 + int(loc[3]) + ((ord(loc[5]) - 65) + x) / 24
    return lat_pos_n, -lon_pos_w  # convert W-positive to E-positive


def latlon_to_locator(lat_deg, lon_deg):
    """Convert lat/lon in decimal degrees (E-positive) to 6-char Maidenhead."""
    if lat_deg == 90:
        return "N POLE"
    if lat_deg == -90:
        return "S POLE"
    lon_w = -lon_deg
    if lon_w == -180:
        lon_w = 180
    zlo = (180 - lon_w) / 20
    zla = (lat_deg + 90) / 10
    za, zb = int(zlo), int(zla)
    zlo, zla = (zlo - za) * 10, (zla - zb) * 10
    zc, zd = int(zlo), int(zla)
    return (chr(65 + za) + chr(65 + zb) + chr(48 + zc) + chr(48 + zd)
            + chr(97 + int((zlo - zc) * 24)) + chr(97 + int((zla - zd) * 24)))


def dms_to_dd(deg, minutes, seconds, direction):
    """(deg, min, sec, 'N'/'S'/'E'/'W') -> signed decimal degrees."""
    dd = abs(deg) + abs(minutes) / 60 + abs(seconds) / 3600
    if direction.upper() in ("S", "W"):
        dd = -dd
    return dd


def dd_to_dms(dd, is_lat=True):
    """Signed decimal degrees -> (deg, min, sec, direction)."""
    if is_lat:
        direction = "N" if dd >= 0 else "S"
    else:
        direction = "E" if dd >= 0 else "W"
    a = abs(dd)
    d = int(a)
    m_full = (a - d) * 60
    m = int(m_full)
    s = int((m_full - m) * 60 + 0.5)
    if s == 60:
        s = 0
        m += 1
    if m == 60:
        m = 0
        d += 1
    return d, m, s, direction


# =============================================================================
# Astronomy: Moon, Sun, sidereal time, az/el
# =============================================================================
def _orbital_constants(datej, t):
    """Compute D1..D19 tracking algorithm constants (VB subroutine 4060)."""
    t5 = datej - 53997.5 + t
    return {
        "t5": t5,
        "d1": fnc(0.606434 + 0.03660110129 * t5),
        "d2": fnc(0.374897 + 0.03629164709 * t5),
        "d3": fnc(0.259091 + 0.03674819520 * t5),
        "d4": fnc(0.827362 + 0.03386319198 * t5),
        "d5": fnc(0.347343 - 0.00014709391 * t5),
        "d7": fnc(0.779072 + 0.00273790931 * t5),
        "d8": fnc(0.993126 + 0.00273777850 * t5),
        "d12": fnc(0.505498 + 0.00445046867 * t5),
        "d13": fnc(0.140023 + 0.00445036173 * t5),
        "d16": fnc(0.053856 + 0.00145561327 * t5),
        "d19": fnc(0.056531 + 0.00023080893 * t5),
    }


def moon_position(datej, t):
    """Geocentric Moon RA, DEC, range (km), and semidiameter (arcmin).

    datej: Julian date offset (BJUL = FNJULIAN(year, month, day))
    t:     day fraction in [0, 1)
    """
    c = _orbital_constants(datej, t)
    t5 = c["t5"]
    d1, d2, d3, d4, d5 = c["d1"], c["d2"], c["d3"], c["d4"], c["d5"]
    d7, d8, d12 = c["d7"], c["d8"], c["d12"]
    sin, cos = math.sin, math.cos

    # Moon ecliptic longitude
    dlon = (22640 * sin(d2) - 4586 * sin(d2 - 2*d4) + 2370 * sin(2*d4)
            + 769 * sin(2*d2) - 668 * sin(d8) - 412 * sin(2*d3)
            - 212 * sin(2*d2 - 2*d4) - 206 * sin(d2 - 2*d4 + d8)
            + 192 * sin(d2 + 2*d4))
    dlon += (165 * sin(2*d4 - d8) + 148 * sin(d2 - d8) - 125 * sin(d4)
             - 110 * sin(d2 + d8) - 55 * sin(2*d3 - 2*d4) - 45 * sin(d2 + 2*d3)
             + 40 * sin(d2 - 2*d3) - 38 * sin(d2 - 4*d4) + 36 * sin(3*d2)
             - 31 * sin(2*d2 - 4*d4) + 28 * sin(d2 - 2*d4 - d8)
             - 24 * sin(2*d4 + d8) + 19 * sin(d2 - d4) + 18 * sin(d4 + d8))
    dlon += (15 * sin(d2 + 2*d4 - d8) + 14 * sin(2*d2 + 2*d4) + 14 * sin(4*d4)
             - 13 * sin(3*d2 - 2*d4))
    dlon += (-11 * sin(d2 + 16*c["d7"] - 18*d12) + 10 * sin(2*d2 - d8)
             + 9 * sin(d2 - 2*d3 - 2*d4) + 9 * cos(d2 + 16*c["d7"] - 18*d12)
             - 9 * sin(2*d2 - 2*d4 + d8) - 8 * sin(d2 + d4)
             + 8 * sin(2*d4 - 2*d8) - 8 * sin(2*d2 + d8))
    dlon += (-7 * sin(2*d8) - 7 * sin(d2 - 2*d4 + 2*d8) + 7 * sin(d5)
             - 6 * sin(d2 - 2*d3 + 2*d4) - 6 * sin(2*d3 + 2*d4)
             - 4 * sin(d2 - 4*d4 + d8)
             + 4 * (t5 / 36525 + 1) * cos(d2 + 16*c["d7"] - 18*d12)
             - 4 * sin(2*d2 + 2*d3))
    dlon += (4 * (t5 / 36525 + 1) * sin(d2 + 16*c["d7"] - 18*d12)
             + 3 * sin(d2 - 3*d4) - 3 * sin(d2 + 2*d4 + d8)
             - 3 * sin(2*d2 - 4*d4 + d8) + 3 * sin(d2 - 2*d8)
             + 3 * sin(d2 - 2*d4 - 2*d8) - 2 * sin(2*d2 - 2*d4 - d8))
    dlon += (-2 * sin(2*d3 - 2*d4 + d8) + 2 * sin(d2 + 4*d4) + 2 * sin(4*d2)
             + 2 * sin(4*d4 - d8) + 2 * sin(2*d2 - d4))
    dlon = d1 + TUPI * dlon / 1296000

    # Moon ecliptic latitude
    dlat = (18461 * sin(d3) + 1010 * sin(d2 + d3) + 1000 * sin(d2 - d3)
            - 624 * sin(d3 - 2*d4) - 199 * sin(d2 - d3 - 2*d4)
            - 167 * sin(d2 + d3 - 2*d4) + 117 * sin(d3 + 2*d4)
            + 62 * sin(2*d2 + d3) + 33 * sin(d2 - d3 + 2*d4)
            + 32 * sin(2*d2 - d3) - 30 * sin(d3 - 2*d4 + d8))
    dlat += (-16 * sin(2*d2 + d3 - 2*d4) + 15 * sin(d2 + d3 + 2*d4)
             + 12 * sin(d3 - 2*d4 - d8) - 9 * sin(d2 - d3 - 2*d4 + d8)
             - 8 * sin(d3 + d5) + 8 * sin(d3 + 2*d4 - d8)
             - 7 * sin(d2 + d3 - 2*d4 + d8) + 7 * sin(d2 + d3 - d8)
             - 7 * sin(d2 + d3 - 4*d4) - 6 * sin(d3 + d8))
    dlat += (-6 * sin(3*d3) + 6 * sin(d2 - d3 - d8) - 5 * sin(d3 + d4)
             - 5 * sin(d2 + d3 + d8) - 5 * sin(d2 - d3 + d8) + 5 * sin(d3 - d8)
             + 5 * sin(d3 - d4) + 4 * sin(3*d2 + d3) - 4 * sin(d3 - 4*d4))
    dlat += (-3 * sin(d2 - d3 - 4*d4) + 3 * sin(d2 - 3*d3)
             - 2 * sin(2*d2 - d3 - 4*d4) - 2 * sin(3*d3 - 2*d4)
             + 2 * sin(2*d2 - d3 + 2*d4) + 2 * sin(d2 - d3 + 2*d4 - d8)
             + 2 * sin(2*d2 - d3 - 2*d4) + 2 * sin(3*d2 - d3))
    dlat = TUPI * dlat / 1296000

    # Convert to equatorial RA/DEC (obliquity ~23.4393 deg: sin=0.397821, cos=0.917463)
    dec1 = cos(dlat) * sin(dlon) * 0.397821 + sin(dlat) * 0.917463
    dec = math.atan2(dec1, math.sqrt(1 - dec1 * dec1))
    rac = cos(dlat) * cos(dlon) / cos(dec)
    ras = (cos(dlat) * sin(dlon) * 0.917463 - sin(dlat) * 0.397821) / cos(dec)
    ra = math.atan2(ras, rac)
    if ra < 0:
        ra += TUPI

    # Earth-Moon distance (in equatorial Earth radii)
    rng = (60.36298 - 3.27746 * cos(d2) - 0.57994 * cos(d2 - 2*d4)
           - 0.46357 * cos(2*d4) - 0.08904 * cos(2*d2)
           + 0.03865 * cos(2*d2 - 2*d4) - 0.03237 * cos(2*d4 - d8)
           - 0.02688 * cos(d2 + 2*d4) - 0.02358 * cos(d2 - 2*d4 + d8)
           - 0.02030 * cos(d2 - d8) + 0.01719 * cos(d4) + 0.01671 * cos(d2 + d8))
    rng += (0.01247 * cos(d2 - 2*d3) + 0.00704 * cos(d8) + 0.00529 * cos(2*d4 + d8)
            - 0.00524 * cos(d2 - 4*d4) + 0.00398 * cos(d2 - 2*d4 - d8)
            - 0.00366 * cos(3*d2) - 0.00295 * cos(2*d2 - 4*d4)
            - 0.00263 * cos(d4 + d8) + 0.00249 * cos(3*d2 - 2*d4)
            - 0.00221 * cos(d2 + 2*d4 - d8))
    rng += (0.00185 * cos(2*d3 - 2*d4) - 0.00161 * cos(2*d4 - 2*d8)
            + 0.00147 * cos(d2 + 2*d3 - 2*d4) - 0.00142 * cos(4*d4)
            + 0.00139 * cos(2*d2 - 2*d4 + d8) - 0.00118 * cos(d2 - 4*d4 + d8)
            - 0.00116 * cos(2*d2 + 2*d4) - 0.00110 * cos(2*d2 - d8))
    semidia = 936.74867 / rng        # arcminutes
    range_km = rng * EARTH_RADIUS_KM
    return ra, dec, range_km, semidia


def sun_position(datej, t):
    """Geocentric Sun RA, DEC at fractional day t of julian date datej."""
    c = _orbital_constants(datej, t)
    t5 = c["t5"]
    d1, d7, d8 = c["d1"], c["d7"], c["d8"]
    d13, d16, d19 = c["d13"], c["d16"], c["d19"]
    sin, cos = math.sin, math.cos

    sunlon = (6910 * sin(d8) + 72 * sin(2*d8) - 17 * (t5/36525 + 1) * sin(d8)
              - 7 * cos(d8 - d19) + 6 * sin(d1 - d7)
              + 5 * sin(4*d8 - 8*d16 + 3*d19) - 5 * cos(2*d8 - 2*d13)
              - 4 * sin(d8 - d13) + 4 * cos(4*d8 - 8*d16 + 3*d19)
              + 3 * sin(2*d8 - 2*d13) - 3 * sin(d19)
              - 3 * sin(2*d8 - 2*d19))
    sunlon = d7 + TUPI * sunlon / 1296000
    sundec_sin = sin(sunlon) * 0.397821
    sundec = math.atan2(sundec_sin, math.sqrt(1 - sundec_sin * sundec_sin))
    sunrac = cos(sunlon) / cos(sundec)
    sunras = sin(sunlon) * 0.917463 / cos(sundec)
    sunra = math.atan2(sunras, sunrac)
    sunra -= int(sunra / TUPI) * TUPI
    if sunra < 0:
        sunra += TUPI
    return sunra, sundec


def gast_hours(datej, t):
    """Greenwich Apparent Sidereal Time in hours [0, 24)."""
    t3 = datej - 35735
    a = 0.0657098232 * t3
    gmst = 6.67170278 + (a - int(a / 24) * 24) + 1.0027379093 * t * 24
    omega = (372.1133 - 0.0529539 * (t3 + t)) * RAD
    sgn = math.copysign(1, omega) if omega != 0 else 1
    omega -= int(sgn * omega / TUPI) * TUPI * sgn
    if abs(omega) < PI:
        omega -= TUPI * sgn
    g = gmst + 0.00029 * math.sin(omega)
    g -= int(g / 24) * 24
    return g


def az_el(ra, dec, gast_hrs, lat_rad, lon_rad_west, range_km=None, parallax=False):
    """Compute geocentric azimuth and elevation (radians) plus GHA (radians).

    lon_rad_west: longitude in radians, positive west (matches original).
    parallax: if True and range_km given, apply geocentric parallax correction.
    """
    gha = gast_hrs * 0.2617994 - ra
    if gha < 0:
        gha += TUPI
    if gha > TUPI:
        gha -= TUPI
    uha = lon_rad_west - gha
    elsin = math.cos(lat_rad) * math.cos(uha) * math.cos(dec) + math.sin(dec) * math.sin(lat_rad)
    elcos = math.sqrt(max(0, 1 - elsin * elsin))
    if elcos == 0:
        return 0.0, 0.0, gha
    el = math.atan2(elsin, elcos)

    az_cos = (math.sin(dec) / (math.cos(lat_rad) * math.cos(el))
              - math.sin(lat_rad) / math.cos(lat_rad) * (math.sin(el) / math.cos(el)))
    az_sin_t = math.sin(lat_rad) * math.sin(dec) + math.cos(lat_rad) * math.cos(dec) * math.cos(uha)
    az_sin = math.sin(uha) * math.cos(dec) / math.sqrt(max(1e-30, 1 - az_sin_t * az_sin_t))
    az = math.atan2(az_sin, az_cos)
    if az <= 0:
        az += TUPI

    if parallax and range_km:
        faz = az - PI
        fel = PI / 2 - el
        phi = math.sin(2 * lat_rad) * 0.000536775 - 0.000000903 * math.sin(4 * lat_rad)
        rho = (0.99832005 + 0.00168349 * math.cos(2 * lat_rad)
               - 0.00000355 * math.cos(4 * lat_rad)
               + 0.00000001 * math.cos(6 * lat_rad))
        h0 = EARTH_RADIUS_KM / range_km
        h0 = math.atan2(math.sin(h0), math.sqrt(1 - math.sin(h0) ** 2))
        az = faz + rho * h0 * math.sin(phi) * math.sin(faz) / math.sin(fel)
        gamma = phi * math.cos(faz)
        el = rho * h0 * math.sin(fel - gamma) + fel
        az += PI
        if az >= TUPI:
            az -= TUPI
        el = PI / 2 - el

    return az, el, gha


def topocentric_moon(lat_rad, ra, gha, dec, range_km, gast_hrs):
    """Compute topocentric Moon RA, GHA, and DEC (radians) for observer at lat_rad."""
    mpar = math.atan2(math.sin(EARTH_RADIUS_KM / range_km),
                      math.cos(math.sqrt(range_km * range_km - EARTH_RADIUS_KM ** 2) / range_km))
    gclat = lat_rad - (0.1924 * TUPI / 360) * math.sin(2 * lat_rad)
    gdist = 0.99833 + 0.00167 * math.cos(2 * lat_rad)
    gaux = math.atan2(math.tan(gclat), math.cos(gha))

    rat = ra - mpar * gdist * math.cos(gclat) * math.sin(gha) / math.cos(dec)
    dect = dec - mpar * gdist * math.sin(gclat) * math.sin(gaux - dec) / math.sin(gaux)
    ghat = gast_hrs / 24 * TUPI - rat
    if ghat < 0:
        ghat += TUPI
    if ghat >= TUPI:
        ghat -= TUPI
    return rat, ghat, dect


# =============================================================================
# Great circle distance between two stations
# =============================================================================
def great_circle_km(lat1, lon1, lat2, lon2):
    """All inputs in radians; returns km along Earth's surface."""
    s = math.sqrt(math.sin((lat1 - lat2) / 2) ** 2
                  + math.cos(lat1) * math.cos(lat2) * math.sin((lon1 - lon2) / 2) ** 2)
    return 2 * math.atan2(s, math.sqrt(1 - s * s)) * DEG * 111.12


# =============================================================================
# Tracking schedule generator
# =============================================================================
def generate_schedule(
    *,
    object_name,
    start_date,
    end_date,
    increment_minutes,
    home,         # dict: lat_dd, lon_dd, locator, callsign, min_el, max_el
    dx=None,      # same shape or None
    units="km",
    region=1,     # 1=Europe/Africa, 2=Americas, 3=Asia/Oceania (for moon windows)
    ra_hours=None, dec_deg=None,  # for non-moon/sun objects
    max_lines=2000,
):
    """Generate a list of tracking schedule rows from start_date to end_date.

    Each row is a dict with timestamp, az/el for each enabled station, polarization,
    moon range, GHA, DEC, etc.
    """
    rows = []
    home_lat = home["lat_dd"] * RAD
    home_lon_w = -home["lon_dd"] * RAD            # original program is W-positive
    home_min = home["min_el"] * RAD
    home_max = home["max_el"] * RAD

    if dx:
        dx_lat = dx["lat_dd"] * RAD
        dx_lon_w = -dx["lon_dd"] * RAD
        dx_min = dx["min_el"] * RAD
        dx_max = dx["max_el"] * RAD
        dx_km = great_circle_km(home_lat, home_lon_w, dx_lat, dx_lon_w)
    else:
        dx_km = None

    one_day = timedelta(days=1)
    cur = start_date
    inc = max(1, int(increment_minutes))

    while cur <= end_date and len(rows) < max_lines:
        datej = fnjulian(cur.year, cur.month, cur.day)

        # Iterate through the day in increments
        total_minutes = 0
        while total_minutes < 1440 and len(rows) < max_lines:
            t = total_minutes / 1440
            gast = gast_hours(datej, t)

            if object_name == "Moon":
                ra, dec, range_km, semidia = moon_position(datej, t)
                parallax = True
            elif object_name == "Sun":
                ra, dec = sun_position(datej, t)
                range_km, semidia = None, None
                parallax = False
            else:
                # User-supplied celestial coordinates (fixed)
                ra = (ra_hours or 0) * TUPI / 24
                dec = (dec_deg or 0) * RAD
                range_km, semidia = None, None
                parallax = False

            haz, hel, gha = az_el(ra, dec, gast, home_lat, home_lon_w, range_km, parallax)
            row = {
                "utc": f"{cur.isoformat()} {total_minutes // 60:02d}:{total_minutes % 60:02d}",
                "home_az": round(fna(haz), 1),
                "home_el": round(fna(hel), 1),
                "gha": round(fna(gha), 1),
                "dec": round(fna(dec), 1),
                "ra_hours": round(ra * 24 / TUPI, 3),
            }

            if range_km is not None:
                rng = range_km if units == "km" else range_km * 0.621371
                row["range"] = round(rng)
                row["semidia_arcmin"] = round(semidia, 2)
                # Topocentric
                rat, ghat, dect = topocentric_moon(home_lat, ra, gha, dec, range_km, gast)
                row["gha_topo"] = round(fna(ghat), 1)
                row["dec_topo"] = round(fna(dect), 1)

            home_visible = home_min <= hel <= home_max

            if dx:
                yaz, yel, _ = az_el(ra, dec, gast, dx_lat, dx_lon_w, range_km, parallax)
                row["dx_az"] = round(fna(yaz), 1)
                row["dx_el"] = round(fna(yel), 1)
                dx_visible = dx_min <= yel <= dx_max
                # Cross polarization angle (only valid when both above horizon)
                if home_visible and dx_visible:
                    pol = fnpol(haz, hel, home_lat) - fnpol(yaz, yel, dx_lat)
                    if abs(pol) >= 90:
                        pol = pol - (180 if pol > 0 else -180)
                    row["pol_deg"] = round(pol)
                visible = home_visible and dx_visible
            else:
                visible = home_visible

            if visible:
                rows.append(row)

            total_minutes += inc
        cur += one_day

    return {
        "rows": rows,
        "dx_path_km": round(dx_km, 1) if dx_km else None,
        "units": units,
        "truncated": len(rows) >= max_lines,
    }


# =============================================================================
# Flask web app
# =============================================================================
app = Flask(__name__)


INDEX_HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>GJTRACKER &mdash; Moon/Sun Tracker</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 1100px; margin: 1em auto; padding: 0 1em; color: #222; }
  h1 { margin-bottom: 0.1em; }
  .sub { color: #666; margin-top: 0; }
  fieldset { border: 1px solid #ccc; border-radius: 4px; margin-bottom: 1em; padding: 0.6em 1em; }
  legend { font-weight: 600; padding: 0 0.5em; }
  .row { display: flex; gap: 1em; flex-wrap: wrap; }
  .col { flex: 1; min-width: 280px; }
  label { display: block; margin-top: 0.4em; font-size: 0.9em; color: #444; }
  input[type=text], input[type=number], input[type=date], select {
    width: 100%; box-sizing: border-box; padding: 0.3em 0.4em;
    border: 1px solid #bbb; border-radius: 3px; font-size: 1em;
  }
  .inline { display: flex; gap: 0.4em; align-items: center; }
  .inline > * { flex: 1; }
  button { padding: 0.5em 1.4em; font-size: 1em; cursor: pointer;
           background: #0066cc; color: white; border: 0; border-radius: 4px; }
  button.secondary { background: #888; }
  button:hover { opacity: 0.9; }
  table { border-collapse: collapse; width: 100%; margin-top: 1em; font-size: 0.85em; }
  th, td { padding: 4px 8px; border-bottom: 1px solid #eee; text-align: right; }
  th { background: #f5f5f5; position: sticky; top: 0; }
  td:first-child, th:first-child { text-align: left; font-family: monospace; }
  .meta { background: #f8f8f8; padding: 0.6em 1em; border-radius: 4px; margin-bottom: 0.5em; }
  .err { color: #c00; }
  .muted { color: #888; font-size: 0.85em; }
</style>
</head>
<body>
<h1>GJTRACKER</h1>
<p class="sub">Moon &amp; Sun position tracker for EME &mdash; Python port of W7GJ's VB6 tool.</p>

<form id="form">
  <fieldset>
    <legend>Object &amp; Time</legend>
    <div class="row">
      <div class="col">
        <label>Object
          <select name="object_name" id="object_name">
            <option>Moon</option>
            <option>Sun</option>
            <option value="Custom">Other (RA/DEC)</option>
          </select>
        </label>
      </div>
      <div class="col" id="custom_radec" style="display:none">
        <label>RA (hours, decimal)
          <input type="number" step="0.001" name="ra_hours" value="0">
        </label>
        <label>DEC (degrees)
          <input type="number" step="0.1" name="dec_deg" value="0">
        </label>
      </div>
      <div class="col">
        <label>Start date <input type="date" name="start_date" required></label>
        <label>End date <input type="date" name="end_date" required></label>
      </div>
      <div class="col">
        <label>Increment (minutes)
          <input type="number" name="increment_minutes" value="15" min="1" max="1440">
        </label>
        <label>Units
          <select name="units"><option>km</option><option>mi</option></select>
        </label>
      </div>
    </div>
  </fieldset>

  <fieldset>
    <legend>Home Station</legend>
    <div class="row">
      <div class="col">
        <label>Callsign <input type="text" name="home_callsign" value="HOME"></label>
        <label>Maidenhead locator
          <input type="text" name="home_locator" placeholder="e.g. CN87xq">
        </label>
        <button type="button" class="secondary" onclick="fromLoc('home')">Locator &rarr; Lat/Lon</button>
        <button type="button" class="secondary" onclick="toLoc('home')">Lat/Lon &rarr; Locator</button>
      </div>
      <div class="col">
        <label>Latitude (decimal degrees, N positive)
          <input type="number" step="0.000001" name="home_lat" value="47.6">
        </label>
        <label>Longitude (decimal degrees, E positive)
          <input type="number" step="0.000001" name="home_lon" value="-122.3">
        </label>
      </div>
      <div class="col">
        <label>Min elevation (deg) <input type="number" name="home_min_el" value="0"></label>
        <label>Max elevation (deg) <input type="number" name="home_max_el" value="90"></label>
      </div>
    </div>
  </fieldset>

  <fieldset>
    <legend><label style="display:inline"><input type="checkbox" id="dx_enabled" name="dx_enabled" checked> DX Station (optional)</label></legend>
    <div class="row" id="dx_panel">
      <div class="col">
        <label>Callsign <input type="text" name="dx_callsign" value="DX"></label>
        <label>Maidenhead locator
          <input type="text" name="dx_locator" placeholder="e.g. JO50">
        </label>
        <button type="button" class="secondary" onclick="fromLoc('dx')">Locator &rarr; Lat/Lon</button>
        <button type="button" class="secondary" onclick="toLoc('dx')">Lat/Lon &rarr; Locator</button>
      </div>
      <div class="col">
        <label>Latitude (decimal degrees, N positive)
          <input type="number" step="0.000001" name="dx_lat" value="50.0">
        </label>
        <label>Longitude (decimal degrees, E positive)
          <input type="number" step="0.000001" name="dx_lon" value="8.0">
        </label>
      </div>
      <div class="col">
        <label>Min elevation (deg) <input type="number" name="dx_min_el" value="0"></label>
        <label>Max elevation (deg) <input type="number" name="dx_max_el" value="90"></label>
      </div>
    </div>
  </fieldset>

  <button type="submit">Calculate Tracking Schedule</button>
  <button type="button" class="secondary" onclick="loadNow()">Use today &rarr; +1 day</button>
</form>

<div id="output"></div>

<script>
const $ = sel => document.querySelector(sel);

// Default dates: today + 1 day
function loadNow() {
  const today = new Date();
  const tomorrow = new Date(today); tomorrow.setDate(today.getDate() + 1);
  $('[name=start_date]').value = today.toISOString().slice(0,10);
  $('[name=end_date]').value = tomorrow.toISOString().slice(0,10);
}
loadNow();

$('#object_name').addEventListener('change', e => {
  $('#custom_radec').style.display = e.target.value === 'Custom' ? '' : 'none';
});

$('#dx_enabled').addEventListener('change', e => {
  $('#dx_panel').style.opacity = e.target.checked ? '1' : '0.4';
});

async function fromLoc(which) {
  const loc = $(`[name=${which}_locator]`).value.trim();
  if (!loc) return alert('Enter a locator first.');
  const r = await fetch('/locator?loc=' + encodeURIComponent(loc));
  const data = await r.json();
  if (data.error) return alert(data.error);
  $(`[name=${which}_lat]`).value = data.lat.toFixed(6);
  $(`[name=${which}_lon]`).value = data.lon.toFixed(6);
}

async function toLoc(which) {
  const lat = parseFloat($(`[name=${which}_lat]`).value);
  const lon = parseFloat($(`[name=${which}_lon]`).value);
  const r = await fetch(`/latlon?lat=${lat}&lon=${lon}`);
  const data = await r.json();
  if (data.error) return alert(data.error);
  $(`[name=${which}_locator]`).value = data.locator;
}

$('#form').addEventListener('submit', async e => {
  e.preventDefault();
  const fd = new FormData(e.target);
  const payload = Object.fromEntries(fd.entries());
  payload.dx_enabled = $('#dx_enabled').checked;
  $('#output').innerHTML = '<p class="muted">Computing&hellip;</p>';

  const r = await fetch('/calculate', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(payload),
  });
  const data = await r.json();
  if (data.error) {
    $('#output').innerHTML = `<p class="err">${data.error}</p>`;
    return;
  }
  renderResults(data);
});

function renderResults(data) {
  if (!data.rows.length) {
    $('#output').innerHTML = '<p class="err">No time slots had the object visible from the selected station(s). Try widening elevation limits or extending the date range.</p>';
    return;
  }
  let html = `<div class="meta">
    Computed <b>${data.rows.length}</b> rows${data.truncated ? ' (truncated)' : ''}.
    ${data.dx_path_km != null ? `Path distance home&hellip;DX: <b>${data.dx_path_km}</b> km` : ''}
  </div>`;

  const r0 = data.rows[0];
  const hasDx = 'dx_az' in r0;
  const hasMoon = 'range' in r0;

  html += '<table><thead><tr>';
  html += '<th>UTC</th><th>Home Az</th><th>Home El</th>';
  if (hasDx) html += '<th>DX Az</th><th>DX El</th><th>Pol &deg;</th>';
  html += '<th>GHA</th><th>DEC</th><th>RA (h)</th>';
  if (hasMoon) html += `<th>Range (${data.units})</th><th>SemiDia &prime;</th><th>GHA topo</th><th>DEC topo</th>`;
  html += '</tr></thead><tbody>';

  for (const row of data.rows) {
    html += `<tr><td>${row.utc}</td><td>${row.home_az}</td><td>${row.home_el}</td>`;
    if (hasDx) html += `<td>${row.dx_az}</td><td>${row.dx_el}</td><td>${row.pol_deg ?? ''}</td>`;
    html += `<td>${row.gha}</td><td>${row.dec}</td><td>${row.ra_hours}</td>`;
    if (hasMoon) html += `<td>${row.range ?? ''}</td><td>${row.semidia_arcmin ?? ''}</td><td>${row.gha_topo ?? ''}</td><td>${row.dec_topo ?? ''}</td>`;
    html += '</tr>';
  }
  html += '</tbody></table>';
  $('#output').innerHTML = html;
}
</script>
</body>
</html>
"""


@app.route("/")
def index():
    return render_template_string(INDEX_HTML)


@app.route("/locator")
def locator_endpoint():
    loc = request.args.get("loc", "").strip()
    try:
        lat, lon = locator_to_latlon(loc)
        return jsonify({"lat": lat, "lon": lon})
    except (ValueError, IndexError) as e:
        return jsonify({"error": f"Invalid locator: {e}"}), 400


@app.route("/latlon")
def latlon_endpoint():
    try:
        lat = float(request.args.get("lat", "0"))
        lon = float(request.args.get("lon", "0"))
        return jsonify({"locator": latlon_to_locator(lat, lon)})
    except ValueError as e:
        return jsonify({"error": str(e)}), 400


@app.route("/calculate", methods=["POST"])
def calculate_endpoint():
    p = request.get_json(force=True)
    try:
        start = date.fromisoformat(p["start_date"])
        end = date.fromisoformat(p["end_date"])
        if end < start:
            return jsonify({"error": "End date must be on or after start date"}), 400
        if (end - start).days > 366:
            return jsonify({"error": "Date range limited to 366 days"}), 400

        home = {
            "lat_dd": float(p["home_lat"]),
            "lon_dd": float(p["home_lon"]),
            "locator": p.get("home_locator", ""),
            "callsign": p.get("home_callsign", ""),
            "min_el": float(p.get("home_min_el", 0)),
            "max_el": float(p.get("home_max_el", 90)),
        }
        dx = None
        if p.get("dx_enabled"):
            dx = {
                "lat_dd": float(p["dx_lat"]),
                "lon_dd": float(p["dx_lon"]),
                "locator": p.get("dx_locator", ""),
                "callsign": p.get("dx_callsign", ""),
                "min_el": float(p.get("dx_min_el", 0)),
                "max_el": float(p.get("dx_max_el", 90)),
            }

        result = generate_schedule(
            object_name=p.get("object_name", "Moon"),
            start_date=start,
            end_date=end,
            increment_minutes=int(p.get("increment_minutes", 15)),
            home=home,
            dx=dx,
            units=p.get("units", "km"),
            ra_hours=float(p.get("ra_hours", 0) or 0),
            dec_deg=float(p.get("dec_deg", 0) or 0),
        )
        return jsonify(result)
    except (KeyError, ValueError) as e:
        return jsonify({"error": str(e)}), 400


if __name__ == "__main__":
    import threading
    import webbrowser

    url = "http://127.0.0.1:5000/"
    threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    print(f"GJTRACKER running at {url}  (Ctrl+C to quit)")
    app.run(host="127.0.0.1", port=5000, debug=False)
