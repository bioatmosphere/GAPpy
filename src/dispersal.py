"""
Inter-site seed dispersal module for GAPpy.

Implements distance-based seed dispersal between neighboring sites using
a negative exponential kernel: dispersal_weight = exp(-distance / max_dispersal_dist).
"""

import math
import numpy as np


# Earth radius in km for haversine calculation
EARTH_RADIUS_KM = 6371.0


def haversine_distance(lat1, lon1, lat2, lon2):
    """Calculate great-circle distance between two points in km.

    Args:
        lat1, lon1: Latitude/longitude of point 1 in decimal degrees.
        lat2, lon2: Latitude/longitude of point 2 in decimal degrees.

    Returns:
        Distance in kilometers.
    """
    lat1_r = math.radians(lat1)
    lat2_r = math.radians(lat2)
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)

    a = (math.sin(dlat / 2.0) ** 2 +
         math.cos(lat1_r) * math.cos(lat2_r) * math.sin(dlon / 2.0) ** 2)
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))

    return EARTH_RADIUS_KM * c


def compute_distance_matrix(sites):
    """Compute pairwise distance matrix between sites.

    Args:
        sites: List of SiteData objects with latitude/longitude attributes.

    Returns:
        2D numpy array of distances in km (nsites x nsites).
    """
    n = len(sites)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = haversine_distance(
                sites[i].latitude, sites[i].longitude,
                sites[j].latitude, sites[j].longitude)
            dist[i, j] = d
            dist[j, i] = d
    return dist


def build_species_index_map(sites):
    """Build a mapping from species unique_id to local index for each site.

    Args:
        sites: List of SiteData objects.

    Returns:
        List of dicts, one per site: {unique_id: local_species_index}.
    """
    maps = []
    for site in sites:
        mapping = {}
        for idx, species in enumerate(site.species):
            mapping[species.unique_id] = idx
        maps.append(mapping)
    return maps


def disperse_seeds(sites, distance_matrix, species_index_maps):
    """Transfer seeds between neighboring sites based on distance.

    For each source site, seeds are exported to neighboring sites using a
    negative exponential dispersal kernel weighted by the species'
    max_dispersal_dist parameter. Seeds are only exported for species that
    have mature trees present (avail_spec > 0) in at least one plot.

    The exported seed quantity per species is:
        export = seed_num * dispersal_weight
    where dispersal_weight = exp(-distance / max_dispersal_dist).

    Seeds are added to the receiving site's plot seedbanks.

    Args:
        sites: List of SiteData objects (with plots already updated for
            the current year).
        distance_matrix: 2D numpy array of inter-site distances in km.
        species_index_maps: List of dicts mapping unique_id to local index.
    """
    nsites = len(sites)
    if nsites < 2:
        return

    # For each source site, determine which species have mature trees
    # (avail_spec averaged across plots > 0 means at least some plots have it)
    for src_idx in range(nsites):
        src_site = sites[src_idx]
        if not src_site.plots or not src_site.species:
            continue

        num_src_species = len(src_site.species)
        src_map = species_index_maps[src_idx]

        # Compute average avail_spec across plots for source site
        avail = np.zeros(num_src_species)
        for plot in src_site.plots:
            avail += plot.avail_spec
        avail /= len(src_site.plots)

        # Export seeds to each destination site
        for dst_idx in range(nsites):
            if dst_idx == src_idx:
                continue

            dist = distance_matrix[src_idx, dst_idx]
            dst_map = species_index_maps[dst_idx]

            for sp_idx in range(num_src_species):
                if avail[sp_idx] <= 0.0:
                    continue

                species = src_site.species[sp_idx]
                if species.max_dispersal_dist <= 0.0:
                    continue

                # Skip if distance exceeds 5x the max dispersal distance
                # (kernel value < 0.007, negligible contribution)
                if dist > 5.0 * species.max_dispersal_dist:
                    continue

                # Negative exponential dispersal kernel
                weight = math.exp(-dist / species.max_dispersal_dist)

                # Seeds exported proportional to seed production and kernel
                seed_export = species.seed_num * avail[sp_idx] * weight

                if seed_export < 1e-6:
                    continue

                # Find matching species in destination site
                dst_sp_idx = dst_map.get(species.unique_id)
                if dst_sp_idx is None:
                    continue

                # Add seeds to destination site's plot seedbanks
                # Distribute evenly across all plots
                per_plot = seed_export / len(sites[dst_idx].plots)
                for plot in sites[dst_idx].plots:
                    plot.seedbank[dst_sp_idx] += per_plot
