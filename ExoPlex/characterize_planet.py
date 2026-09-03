#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:24:18 2026

@author: joesch
"""

import requests
import numpy as np
import scipy.stats as sp
import pandas as pd

class abundances:
    
    def get_solar_norm(solar_norm_name = 'lodders09', verbose = False):
        
        # This function retreives A(X) values for each element X for a desired
        # solar normalization. To this end, it queries the Hypatia catalog,
        # which contains the following solar normalizations: 
        # 'asplund05', 'lodders09', 'anders89', 'grevesse98', 
        # 'asplund09', 'grevesse07', 'absolute', 'original'.
        # See Hypatia API page for further details.
        
        get_normalizations = requests.get("https://hypatiacatalog.com/hypatia/api/v2/solarnorm")
        normalizations = get_normalizations.json()
        solar_norm_list = np.array([normalizations[i]['id'] for i in range(0, len(normalizations))])
        solar_norm = normalizations[np.where(solar_norm_list == solar_norm_name)[0][0]]
        norm = {'normalization_name':solar_norm['id']}
        norm.update(solar_norm['values'])
        
        if verbose:
            print(f'Using {norm['normalization_name']} solar normalization:')
            print(f'A(Fe) = {norm['Fe']}')
            print(f'A(Mg) = {norm['Mg']}')
            print(f'A(Si) = {norm['Si']}')

        return norm
    
    def get_single_star_abundances(star_name, solar_norm_name = 'lodders09', verbose = False):
        
        # Queries the Hypatia Catalog for [Fe/H], [Mg/H], and [Si/H] for an input star
        # Returns a dictionary for the star containing:
        #   'in_hypatia' -- True if star is in Hypatia, False if not
        #   FeH: median [Fe/H] abundance in dex
        #   sig_FeH: Hypatia-listed error for [Fe/H] 
        #   MgH: median [Mg/H] abundance in dex
        #   sig_MgH: Hypatia-listed error for [Mg/H] 
        #   SiH: median [Si/H] abundance in dex
        #   sig_SiH: Hypatia-listed error for [Si/H] 
        
        # If any of the abundance parameters for Fe, Mg, Si are not present for
        # the desired star, a NaN is returned for the given parameter.

        params = {"name": star_name, "element": ["fe"], "solarnorm": [solar_norm_name]}
        star_entry_fe = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
        star_entry_fe = star_entry_fe.json()[0]
        star = {'name': star_name,
                'solar_norm': solar_norm_name,
                'FeH': np.nan,
                'sig_FeH': np.nan,
                'MgH': np.nan,
                'sig_MgH': np.nan,
                'SiH': np.nan,
                'sig_SiH': np.nan}
        
        if star_entry_fe['name'] == 'not-found' and verbose == True:
            print(f'{star_name} not found in Hypatia.')     
        else: 
            star['FeH'] = star_entry_fe['median_value']
            star['sig_FeH'] = star_entry_fe['plusminus']
            
            params = {"name": star_name, "element": ["mg"], "solarnorm": [solar_norm_name]}
            star_entry_mg = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
            star_entry_mg = star_entry_mg.json()[0]
            star['MgH'] = star_entry_mg['median_value']
            
            if star['MgH'] != None:    
                star['sig_MgH'] = star_entry_mg['plusminus']
            else:
                star['MgH'] = np.nan

            
            params = {"name": star_name, "element": ["si"], "solarnorm": [solar_norm_name]}
            star_entry_si = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", params=params)
            star_entry_si = star_entry_si.json()[0]
    
            star['SiH'] = star_entry_si['median_value']
            if star['SiH'] != None:    
                star['sig_SiH'] = star_entry_mg['plusminus']
            else:
                star['SiH'] = np.nan
                
            if verbose == True:
                print(f'Hypatia abundance values for {star['name']} relative to {star['solar_norm']}:')
                print(f"[Fe/H] = {star['FeH']:.2f} +/- {star['sig_FeH']:.2f}") 
                print(f"[Mg/H] = {star['MgH']:.2f} +/- {star['sig_MgH']:.2f}") 
                print(f"[Si/H] = {star['SiH']:.2f} +/- {star['sig_SiH']:.2f}") 
        
        return star
    
    def generate_abundance_samples(star, num_samples = 500, verbose = False):
        
        # Generates X number of random [Fe/H], [Mg/H], and [Si/H] samples.
        #
        # If the desired star has [Fe/H], [Mg/H], and [Si/H] and corresponding
        # uncertainties, then the X abundance samples are sample directly from
        # these values assuming they are normally distributed (e.g., 
        # FeH_samples = Normal([Fe/H], sig_[Fe/H], X).
        # 
        # IMPORTANT: If the specified star is listed as galactic_nrpz OR it is 
        # missing a measurement for one or more of the major rock-building elements
        # the [Fe/H], [Mg/H], and [Si/H] values for X random Hypatia stars are 
        # returned instead.
        #
        # Whether the samples are generated directly from a desired host's measured
        # abundances or are a sample of X random Hypatia stars is noted
        # under "Abundance sample provenance" and added to the star dictionary.

        
        if star['name'] != 'galactic_nrpz' and np.isnan(star['MgH']) ==  False and np.isnan(star['FeH']) ==  False:
            if verbose:
                print(f'Generating {num_samples} abundance samples for {star}')
            star.update({'Abundance sample provenance':f'{star['name']} listed abundance values'})
            FeH_samples = sp.norm.rvs(star['FeH'], star['sig_FeH'], num_samples)
            MgH_samples = sp.norm.rvs(star['MgH'], star['sig_MgH'], num_samples)
            SiH_samples = sp.norm.rvs(star['SiH'], star['sig_SiH'], num_samples)
            star.update({'FeH_samples':FeH_samples, 
                         'MgH_samples':MgH_samples,
                         'SiH_samples':SiH_samples})
        else:
            if verbose:
                print(f'Grabbing values for {num_samples} random Hypatia stars with [Fe/H], [Mg/H], AND [Si/H].')
            star.update({'Abundance sample provenance':f'generated from {num_samples} random Hypatia stars'})
            hypatia_stars = pd.read_csv('./Hypatia/hypatia_abundance_sample.csv', comment = '#', index_col = False)
            inds = np.random.randint(0, len(hypatia_stars), num_samples)
            FeH_samples = np.array(hypatia_stars['FeH'])[inds]
            MgH_samples = np.array(hypatia_stars['MgH'])[inds]
            SiH_samples = np.array(hypatia_stars['SiH'])[inds]
            star.update({'FeH_samples':FeH_samples, 
                         'MgH_samples':MgH_samples,
                         'SiH_samples':SiH_samples})
        return star
            
            
            


            
            
        
        
            
        