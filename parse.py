import math
import re
import geopandas as gpd
import pandas as pd
from pyproj import Transformer, CRS
from shapely.geometry import Point, Polygon

excel_dataframe = pd.read_excel('D:\\study\\NextGIS\\source.xlsx')
geopackage_geodataframe = gpd.GeoDataFrame(columns=excel_dataframe.columns.tolist() + ['geometry'])

KEYWORDS = {
    'Окружность (с заданным радиусом)': 'радиус',
    'Пересечение полигонов': 'захоронение',
    'Окружность (с заданными координатами': 'координаты окружности',
    'Вычитание полигонов': 'исключенного',
    'Объединение полигонов': 'используемого',
}

RE_EXP = {
    'pattern_radius':
        r"(\d+(?:\.\d+)?) (\w+)",
    'pattern_coordinates':
        r"ш=(\d+)[^.](\d+(?:\.\d+)?)[^.](\d+(?:\.\d+)?)*\D+(\d+)[^.](\d+(?:\.\d+)?)[^.](\d+(?:\.\d+)?)*[^.]",
    'pattern_coordinates_without_numbers':
        r"(\d+)[^.](\d+(?:\.\d+)?)[^.](\d+(?:\.\d+)?)*\D+(\d+)[^.](\d+(?:\.\d+)?)[^.](\d+(?:\.\d+)?)*"
}

CRS_MAP = {
    "wgs-84": 'EPSG:4326',
    "ск-42": 'EPSG:4284'
}


# Points must have Polygon or Point data type
def crs_transform(points, crs):
    global transformed_points
    target_crs = CRS_MAP["wgs-84"]

    if crs != target_crs:
        source_crs = CRS.from_string(crs)
        target_crs = CRS.from_string(target_crs)
        transformer = Transformer.from_crs(source_crs,
                                           target_crs,
                                           always_xy=True)
        if isinstance(points, Polygon):
            x_coord, y_coord = points.exterior.xy
            x_coords_transformed, y_coords_transformed = transformer.transform(x_coord, y_coord)
            transformed_points = Polygon(zip(x_coords_transformed, y_coords_transformed))
        elif isinstance(points, Point):
            lon, lat = point.x, point.y
            lon_transformed, lat_transformed = transformer.transform(lon, lat)
            transformed_points = Point(lon_transformed, lat_transformed)
        return transformed_points

    else:
        return points


def determinate_type_of_stroke(coordinates):
    if "ш=" in coordinates:
        pattern = RE_EXP['pattern_coordinates']
    else:
        pattern = RE_EXP['pattern_coordinates_without_numbers']

    return pattern


def decimal_coordinates(deg, minute, sec):
    if sec == '':
        sec = 0.0

    if float(deg) > 180:
        deg_int = str(round(float(deg)))
        if float(deg) != 0 and float(minute) != 0 and float(sec) != 0:
            deg = (float(deg) // 10)
        else:
            sec = minute
            if len(deg_int) == 4:
                deg, minute = deg[:2], deg[2:]
            elif len(deg_int) == 5:
                deg, minute = deg[:2], deg[3:]

    if float(minute) > 60 or len(minute) >= 4:
        if float(deg) != 0 and float(minute) != 0 and float(sec) != 0:
            minute = minute[2:]
        else:
            minute_int = str(round(float(minute)))

            if len(minute_int) == 4:
                minute, sec = minute[:2], minute[2:]
            elif len(minute_int) == 5:
                minute, sec = minute[:2], minute[3:]

    return float(deg) + (float(minute) / 60) + (float(sec) / 3600)


def separation_coordinates(stroke, separator):
    stroke = stroke.split(separator)
    return stroke[0], stroke[1]


def creating_polygon(pattern, coordinates, source_crs):
    global polygon

    points_dec = []
    coordinates = re.sub(r"\d{2,}[0]\d{2,}\W\d",
                         lambda match: match.group().replace('0', 'о'),
                         coordinates)
    points = re.findall(pattern, coordinates)

    for p in points:
        lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec = p[:6]
        lat = decimal_coordinates(lat_deg, lat_min, lat_sec)
        lon = decimal_coordinates(lon_deg, lon_min, lon_sec)
        points_dec.append([lon, lat])

    if len(points_dec) >= 3:
        first_point = points_dec[0]
        points_dec.append(first_point)
        polygon = Polygon(points_dec)
    elif len(points_dec) == 1:
        polygon = Point(points_dec)

    crs = source_crs
    polygon = crs_transform(polygon, crs)

    return polygon


def creating_point_or_circle(pattern, coordinates, source_crs, radius=0):
    global point, lat

    coordinates = re.sub(r"\d{2,}[0]\d{2,}\W\d",
                         lambda match: match.group().replace('0', 'о'),
                         coordinates)
    points = re.findall(pattern, coordinates)

    if points:
        lat_deg, lat_min, lat_sec, lon_deg, lon_min, lon_sec = points[0][:6]
        lat = decimal_coordinates(lat_deg, lat_min, lat_sec)
        lon = decimal_coordinates(lon_deg, lon_min, lon_sec)
        point = Point(lon, lat)

    if radius != 0:
        coef = 2 * math.pi * 6371000 * math.cos(math.radians(lat))
        radius_deg = radius / (coef * 360)
        point = point.buffer(radius_deg)

    crs = source_crs
    point = crs_transform(point, crs)

    return point


polygon_burial = []

for i, row in excel_dataframe.iterrows():

    coords_curr = row["Район захоронения донного грунта (Географические координаты)"]

    coords_curr = coords_curr.casefold(). \
        replace('\n', ''). \
        replace('``', "″"). \
        replace("’", "'"). \
        replace("’’", "″"). \
        replace('..', '.'). \
        replace(',', '.'). \
        replace("центр", "радиус"). \
        replace("pulkovo1942", "ск-42"). \
        replace("pulkovo 1942", "ск-42"). \
        replace("пулково ск-42", "ск-42"). \
        replace("пулково 1942", "ск-42"). \
        replace("(", ""). \
        replace(")", "")

    coords_curr = re.sub(r"(\d)(ш)", r"\1.\2", coords_curr)

    # processing CRS information
    if "ск-42" in coords_curr and "wgs-84" not in coords_curr:
        source_crs = CRS_MAP["ск-42"]
    else:
        source_crs = CRS_MAP["wgs-84"]

    if "ск-42" in coords_curr and "wgs-84" in coords_curr:
        if coords_curr.find("ск-42") < coords_curr.find("wgs-84"):
            coords_curr = coords_curr[coords_curr.find("wgs-84"):]
        else:
            coords_curr = coords_curr[:coords_curr.find("ск-42")]
    coords_curr = coords_curr. \
        replace("ск-42", ""). \
        replace("wgs-84", "")

    # geometry information processing
    keyword_circle_with_radius = KEYWORDS['Окружность (с заданным радиусом)']

    keyword_intersection = KEYWORDS['Пересечение полигонов']
    keyword_circle = KEYWORDS['Окружность (с заданными координатами']
    keyword_difference = KEYWORDS['Вычитание полигонов']
    keyword_union = KEYWORDS['Объединение полигонов']

    if keyword_circle_with_radius in coords_curr:
        # determinate_radius
        pattern_radius = RE_EXP['pattern_radius']
        radius_match = re.search(pattern_radius, coords_curr)
        if radius_match:
            if "морск" in radius_match.group(2):
                radius = float(radius_match.group(1)) * 1852
            elif "кбт" in radius_match.group(2):
                radius = float(radius_match.group(1)) * 185.2
            elif "км" in radius_match.group(2):
                radius = float(radius_match.group(1)) * 1000
            elif "м" in radius_match.group(2):
                radius = float(radius_match.group(1))

        coords_curr = coords_curr.replace(' ', '')

        if keyword_intersection in coords_curr:
            coords_curr_center, coords_curr_burial = separation_coordinates(coords_curr, keyword_intersection)
            pattern = determinate_type_of_stroke(coords_curr)
            circle_burial = creating_point_or_circle(pattern, coords_curr_center, source_crs, radius)
            polygon_burial = creating_polygon(pattern, coords_curr_burial, source_crs)
            try:
                if not polygon_burial.is_valid:
                    polygon_burial = polygon_burial.buffer(0)
                if not circle_burial.is_valid:
                    circle_burial = circle_burial.buffer(0)
                polygon_burial = polygon_burial.intersection(circle_burial)
            except Exception as e:
                print("Ошибка при обрезании полигона:", e)


        elif keyword_circle in coords_curr:
            coords_curr_center, coords_curr_burial = separation_coordinates(coords_curr, keyword_circle)
            pattern = determinate_type_of_stroke(coords_curr)
            circle_burial = creating_point_or_circle(pattern, coords_curr_center, source_crs).buffer(0)
            polygon_burial = creating_polygon(pattern, coords_curr_burial, source_crs)
            polygon_burial = polygon_burial.union(circle_burial) if polygon_burial and circle_burial else None

        else:
            pattern = determinate_type_of_stroke(coords_curr)
            polygon_burial = creating_point_or_circle(pattern, coords_curr, source_crs, radius)

    else:

        coords_curr = coords_curr. \
            replace(' ', ''). \
            replace("непроизводитьзахоронение", "исключенного")

        if keyword_difference in coords_curr:
            coords_curr_burial, coords_curr_not_burial = separation_coordinates(coords_curr, keyword_difference)
            pattern = determinate_type_of_stroke(coords_curr)
            polygon_burial = creating_polygon(pattern, coords_curr_burial, source_crs)
            polygon_not_burail = creating_polygon(pattern, coords_curr_not_burial, source_crs)
            polygon_burial = polygon_burial.difference(polygon_not_burail) \
                if polygon_burial and polygon_not_burail else None


        elif keyword_union in coords_curr:
            coords_curr_first, coords_curr_second = separation_coordinates(coords_curr, keyword_union)
            pattern = determinate_type_of_stroke(coords_curr)
            polygon_first = creating_polygon(pattern, coords_curr_first, source_crs)
            polygon_second = creating_polygon(pattern, coords_curr_second, source_crs)
            polygon_burial = polygon_first.union(polygon_second) if polygon_first and polygon_second else None


        else:
            pattern = determinate_type_of_stroke(coords_curr)
            polygon_burial = creating_polygon(pattern, coords_curr, source_crs)

    if polygon_burial:
        geopackage_geodataframe.loc[i] = row.values.tolist() + [polygon_burial]

geopackage_geodataframe.to_file('D:\\study\\NextGIS\\file.gpkg', driver='GPKG', dtype={'geometry': 'Geometry'})
