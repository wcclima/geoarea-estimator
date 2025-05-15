
import numpy as np
from scipy.special import erfinv
from prettytable import PrettyTable

__all__ = ["GeoRegionArea"]

class GeoRegionArea(object):

    """
    Estimates the area of a given region on Earth.

    GeoRegionArea computes the area within a given boundary by sampling 
    points on the surface of a sphere.  

    Attributes:
        n_samples (int):
            The number of samples used to estimate the area.

        region_area (float):
            The estimate value of the area.

        error (float):
            The confidence interval @ 95% confidence level.

            
    Methods:
        estimate(data, n_samples):    
            Estimates the area of a region on Earth within a given boundary using 
            Monte Carlo sampling.

        results():
                Generate a summary table of the final area estimation results.


    Notes:

        - Since Earth is better approximated by a spheroid, we calculate the area of the desired 
        region by locally approximating Earth by a sphere and using the local geocentric radius 
        distance.

        - The geocentric radius distance is calculated from the average minimum and maximum 
        latitudes and longitudes of the given reagion's boundary.

        - The geocentric radius distance is computed considering a spheroid with equatorial radius 
        a =  6378.1370km and polar radius b = 6356.7523km.
        
        - It gives better approximations for the area when area < (geocentric radius)**2.
    
    """

    def __init__(self):
        """
        Initialises the GeoRegionArea class.

        """

        self.n_samples = None
        self.region_area = None
        self.error = None
    

    def _area_square(self, min_lat, max_lat, min_lon, max_lon):
        """
        Estimate the area of a spherical lat/lon-bounded region on Earth's surface using geocentric radius
        as the radius of the sphere.

        Keyword arguments:
        
        min_lat (float):
            Minimum latitude in radians.
        max_lat (float):
            Maximum latitude in radians.
        min_lon (float):
            Minimum longitude in radians.
        max_lon (float):
            Maximum longitude in radians.

        Returns

        area_square : float
            Estimated surface area of the quadrilateral region in square kilometers.

        Notes
        
        - The area is computed by integrating over a spherical cap adjusted for latitude.
        - Inputs must be in radians.
        """

        a_radius_km = 6378.1370
        b_radius_km = 6356.7523
        mean_lat = 0.5*(max_lat + min_lat)

        geocent_radius_km = np.sqrt(
            (a_radius_km**4*np.cos(mean_lat)**2 + b_radius_km**4*np.sin(mean_lat)**2
             )/(a_radius_km**2*np.cos(mean_lat)**2 + b_radius_km**2*np.sin(mean_lat)**2)
             )

        area_square = geocent_radius_km**2*(np.sin(max_lat) - np.sin(min_lat))*(max_lon - min_lon)

        return area_square


    def _3d_vector(self, lat, lon):
        """
        Convert geographic coordinates to a 3D Cartesian unit vector.

        Keyword arguments
        
        lat (float):
            Latitude of the point in degrees.
        lon (float):
            Longitude of the point in degrees.

        Returns
        
        x (float):
            X-component of the unit vector.
        y (float):
            Y-component of the unit vector.
        z (float):
            Z-component of the unit vector.

        Notes
        -----
        - The conversion assumes a spherical Earth.
        - The returned vector lies on the unit sphere (radius = 1).
        - Input angles are converted internally from degrees to radians.
        """

        lat_rad = np.radians(lat)
        lon_rad = np.radians(lon)

        x = np.sin(np.pi/2. - lat_rad)*np.cos(lon_rad)
        y = np.sin(np.pi/2. - lat_rad)*np.sin(lon_rad)
        z = np.cos(np.pi/2. - lat_rad)

        return x, y, z
    

    def _intersect(
            self, 
            vertices_lat, 
            vertices_lon,  
            point_lat, 
            point_lon
            ):
        
        """
        Checks whether the west-east great circle starting at the point p = (point_lat, point_lon)
        crosses the sides of the polygon (i.e. any great circle segment between two consecutive 
        polygon vertices). The number of crossings is determined by checking 
        point_lon with respect to the longitude of the crossing point, if it exists.

        Returns an array-like of bools, where True means a crossing.


        Keyword arguments:
            vertices_{lat/lon} (float):
                The {latitude/longitude} of the polygon vertices in degrees.

            p_{lat/lon} (float):
                The {latitude/longitude} of the point of interest p in degrees.


        Returns: 
            crossed (np.ndarray): ndarray of bools with shape (n_vertices,) 
                Whether the great circle starting at p crossed the sides of the polygon.

        """
        n_vertices = len(vertices_lat) - 1
        vertices_x, vertices_y, vertices_z = self._3d_vector(vertices_lat, vertices_lon)
        p_x, p_y, p_z = self._3d_vector(np.repeat(point_lat, n_vertices), np.repeat(point_lon, n_vertices))
        cos_angle_v1v2 = vertices_x[:n_vertices]*vertices_x[1:] + vertices_y[:n_vertices]*vertices_y[1:] + vertices_z[:n_vertices]*vertices_z[1:]

        epsilon = 1e-15
        
        norm_v1v2 = np.sqrt(1. - cos_angle_v1v2**2)
        norm_v1v2 += (norm_v1v2 < np.finfo(float).eps)*epsilon
        t_12_x = (vertices_x[1:] - cos_angle_v1v2*vertices_x[:n_vertices])/norm_v1v2
        t_12_y = (vertices_y[1:] - cos_angle_v1v2*vertices_y[:n_vertices])/norm_v1v2
        t_12_z = (vertices_z[1:] - cos_angle_v1v2*vertices_z[:n_vertices])/norm_v1v2

        t_p_x = -np.sin(np.pi/2. - np.radians(point_lat))*np.sin(np.radians(point_lon))
        t_p_y = np.sin(np.pi/2. - np.radians(point_lat))*np.cos(np.radians(point_lon))
        t_p_z = 0

        cos_angle_pv1 = p_x*vertices_x[:n_vertices] + p_y*vertices_y[:n_vertices] + p_z*vertices_z[:n_vertices]
        cos_angle_tpv1 = t_p_x*vertices_x[:n_vertices] + t_p_y*vertices_y[:n_vertices] + t_p_z*vertices_z[:n_vertices]
        cos_angle_pt12 = p_x*t_12_x + p_y*t_12_y + p_z*t_12_z
        cos_angle_tpt12 = t_p_x*t_12_x + t_p_y*t_12_y + t_p_z*t_12_z

        alpha_prime = np.arctan2(1. - cos_angle_pv1**2 - cos_angle_pt12**2,cos_angle_tpv1*cos_angle_pv1 + cos_angle_tpt12*cos_angle_pt12)

        cond1 = (vertices_lat[:n_vertices] > point_lat) != (vertices_lat[1:] > point_lat)
        cond2 = np.radians(point_lon) < alpha_prime

        crossed = cond1 & cond2
        
        return crossed       


    def _internal_point(self, vertices_lat, vertices_lon, points_lat, points_lon):

        """
        Checks whether points p are inside the polygon.

        For a given point, it counts how many times the great circle 
        starting at it crosses the sides of the polygon. An odd number of crosses 
        indicates that the point is inside the polygon and an even number of crosses 
        indicates it is outside.

        Returns True if the point is inside the polygon and False otherwise.

        Keyword arguments:
            p (array): array-like of shape (2,) in the case of one point or
            of shape (n_points, 2) in the case of n_points.
                The latitude/longitude of the point or a list of pairs for a 
                number of points

            labels(array, default = None): array-like of shape (n_points,)
                The labels of the points.

        Returns:
            Table of the point labels, their coordinates and whether they 
            are inside or ouside the polygon.
            
        """

        inside = []

        for p_lat, p_lon in zip(points_lat, points_lon):

            crossings = self._intersect(
                vertices_lat, 
                vertices_lon,
                p_lat, 
                p_lon
                ).sum()
            
            if crossings % 2:
                inside.append(True)
            else:
                inside.append(False)

        return np.array(inside)
    
    def estimate(self, data, n_samples = 1000):
        """
        Estimate the area of a region on Earth using the Monte Carlo method.

        This method performs Monte Carlo integration to approximate the area 
        enclosed by a boundary defined by geographic coordinates (latitude and longitude). 

        Keyword arguments:
        
        data (array): array-like of shape (n_points, 2)
            A list or array of [latitude, longitude] pairs (in degrees) defining the 
            closed boundary of the region to be measured.

        n_samples (int, default=1000):
            The number of random points to sample within the bounding box of the region 
            for the Monte Carlo integration.

        Notes
        
        - The method assumes the Earth is an oblate spheroid for computing bounding area.
        - Latitude and longitude values are expected in degrees, in decimals.
        """

        self.n_samples = n_samples

        region_lat = np.array([np.radians(i[0]) for i in data])
        region_lon = np.array([np.radians(i[1]) for i in data])

        min_lat, max_lat = region_lat.min(), region_lat.max()
        min_lon, max_lon = region_lon.min(), region_lon.max()

        samples_lat = (max_lat - min_lat)*np.random.rand(n_samples) + min_lat
        samples_lon = (max_lon - min_lon)*np.random.rand(n_samples) + min_lon

        p_region = self._internal_point(region_lat, region_lon, samples_lat, samples_lon).sum()/n_samples

        a_tot = self._area_square(min_lat, max_lat, min_lon, max_lon)
        self.region_area = p_region*a_tot
        self.error = erfinv(0.95)*a_tot/np.sqrt(2.*n_samples)
    

    def results(self):
        """
        Generate a summary table of the final area estimation results.

        This method compiles and returns a PrettyTable summarizing key information from the 
        Monte Carlo area estimation. It includes the number of samples, the most recent 
        estimate of the area, and the corresponding 95% confidence interval.

        Returns
        -------
        PrettyTable
            A formatted summary table containing:
                - Number of samples used
                - The estimated area
                - 95% confidence interval around the estimate

        Notes
        -----
        This method assumes that the `estimate` method has already been called.
        """

        error_message = "Estimation has not been performed. Please call `estimate` before `results`."
        if self.region_area is None or self.error is None:
            raise ValueError(error_message)

        estimation_table = PrettyTable(["quantity", "value"])
        estimation_table.add_row(["samples", self.n_samples])
        estimation_table.add_row(["area (km^2)", np.round(self.region_area, 4)])
        estimation_table.add_row(["conf. interval @95%", f'[{np.round(self.region_area - self.error, 4)}, {np.round(self.region_area + self.error, 4)}]'])

        return estimation_table