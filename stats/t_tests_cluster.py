#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 10:19:59 2021

@author: leupinv
"""
# functions that need to be revised to run TFCE and cluster analyses
from mne.channels import find_ch_adjacency
from mne.stats import spatio_temporal_cluster_1samp_test,f_mway_rm
from mne.stats import permutation_cluster_1samp_test, f_threshold_mway_rm, spatio_temporal_cluster_test
import mne
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def evoked_plot(evok_effect,tmin):
    biosemi_montage = mne.channels.make_standard_montage('biosemi128')

    info = mne.create_info(ch_names=biosemi_montage.ch_names, sfreq=256.,
                                ch_types='eeg')
    evok=mne.EvokedArray(evok_effect,info,tmin=tmin)
    evok.set_montage(biosemi_montage)
    return evok

def get_DF_TFCE(evoked,crop_t=False,crop_value=None):
    if crop_t==True:
        data_crop=evoked[0].crop(crop_value[0],crop_value[1])
        data_shape=data_crop.data.shape
        subj_len=len(evoked)
        print(data_shape)
    else:
        data_shape=evoked[0].data.shape
        subj_len=len(evoked)

    X=np.empty((subj_len,data_shape[1],data_shape[0]))
    for idx, ev in enumerate(evoked):
        X[idx,:,:]=ev.crop(crop_value[0],crop_value[1]).data.T
    print(X.shape)
    return X


def get_Ttest_TFCE(X,plot_times='peaks',adjacency=None,averages=None,permutations=None):
    tfce = dict(start=0, step=.5)

    t_obs, clusters, cluster_pv, h0 = spatio_temporal_cluster_1samp_test(
    X, tfce, adjacency=adjacency,
    n_permutations=permutations,out_type='mask',n_jobs=-1)  # a more standard number would be 1000+

    significant_points = cluster_pv.reshape(t_obs.shape).T < .05
    print(str(significant_points.sum()) + " points selected by TFCE ...")
    biosemi_montage = mne.channels.make_standard_montage('biosemi128')
    n_channels = len(biosemi_montage.ch_names)
    info = mne.create_info(ch_names=biosemi_montage.ch_names, sfreq=256.,
                                ch_types='eeg')
    evok=mne.EvokedArray(t_obs.T,info,tmin=-0.3)
    evok.set_montage(biosemi_montage)
    evok.plot_image(mask=significant_points,scalings=1,units='T-value',show_names='auto')

    evok.plot_topomap(plot_times,outlines='head',scalings=1,units='T-value',average=averages,mask=significant_points)
    return t_obs, clusters, cluster_pv, h0

def tTest_TFCE_ana(evoked,crop_t=False,crop_value=None,plot_times='peaks',averages=None,permutations=1000,TFCE=False,p_val=0.05):
    adjacency, _ = find_ch_adjacency(evoked[0][0].info, "eeg")
    evoked_1=get_DF_TFCE(evoked[0],crop_t=True,crop_value=crop_value)
    evoked_2=get_DF_TFCE(evoked[1],crop_t=True,crop_value=crop_value)
    X=evoked_1-evoked_2
    if TFCE:
        t_obs, clusters, cluster_pv, h0=get_Ttest_TFCE(X,plot_times=plot_times,averages=averages,permutations=permutations,adjacency=adjacency)
        return t_obs, clusters, cluster_pv, h0
    else:
        t_obs, clusters, cluster_pv, h0=get_Ttest_cluster(X,
        plot_times=plot_times,averages=averages,permutations=permutations,
        adjacency=adjacency,p_val=p_val)
        return t_obs, clusters, cluster_pv, h0

def get_Ttest_cluster(X,plot_times='peaks',adjacency=None,averages=None,permutations=1000,p_val=0.05):
    n_suj = np.shape(X)[0]
    t_threshold = -stats.distributions.t.ppf(p_val / 2., n_suj - 1)
    t_obs, clusters, cluster_pv, h0 = spatio_temporal_cluster_1samp_test(
    X, adjacency=adjacency,
    n_permutations=permutations,out_type='mask', threshold=t_threshold,n_jobs=-1)  # a more standard number would be 1000+

    T_obs_plot=0*np.ones_like(t_obs)
    for c,p_val in zip(clusters,cluster_pv):
        if p_val<=.05:
            T_obs_plot[c]=t_obs[c]
    data=T_obs_plot.T

    #data=np.flipud(data)

    biosemi_montage = mne.channels.make_standard_montage('biosemi128')
    n_channels = len(biosemi_montage.ch_names)
    info = mne.create_info(ch_names=biosemi_montage.ch_names, sfreq=256.,
                                ch_types='eeg')
    evok=mne.EvokedArray(data,info,tmin=-0.1)
    evok.set_montage(biosemi_montage)
    evok.plot_image(scalings=1,units='T-value',show_names='auto')

    evok.plot_topomap(plot_times,outlines='head',scalings=1,units='T-value',average=averages)





    return t_obs, clusters, cluster_pv,T_obs_plot







def clus_Anovas_ana(evoked,effect_label, crop_value=None,g_excl=None,factor_levels=[2,2],
               effects='A:B',FDR=False,report=None,topo_times='peaks',p_val=0.05,plot_average=0.02,png=None,n_perm=1000,TFCE=None):

    '''
    Note: each effect has to be computed separately:
        effects= "A", first main effect
        effects= "B", second main effect
        effects= "A:B", interaction effect

    '''
    adjacency, _ = find_ch_adjacency(evoked[0][0].info, "eeg")
    #get raw X
    X=[get_DF_TFCE(X,crop_t=True,crop_value=crop_value) for X in evoked]
    #format X
    ## define dimensions
    n_rep=len(evoked[0])
    n_conditions=len(evoked)
    n_chan=128
    n_times=evoked[0][0].data.shape[1]

    ## reformat data
    data = np.swapaxes(np.asarray(X), 1, 0)
    # reshape last two dimensions in one mass-univariate observation-vector
    data = data.reshape(n_rep, n_conditions, n_chan * n_times)
    if TFCE==None:
        TFCE=f_threshold_mway_rm(n_rep,factor_levels,effects,pvalue=p_val)

    print(data.shape)

    def stat_fun_anov(*args):
        return f_mway_rm(np.swapaxes(args,1,0),factor_levels=factor_levels,effects=effects,return_pvals=False)[0]

    # Compute Anova
    F_obs,clusters,clust_p,H0=clu=\
    spatio_temporal_cluster_test(X,adjacency=adjacency,threshold=TFCE,stat_fun=stat_fun_anov,n_permutations=n_perm)
    #fvals, pvals = f_mway_rm(data, factor_levels, effects=effects,correction=True)

    #Plot anova

    F_obs_plot=0*np.ones_like(F_obs)
    for c,p_valu in zip(clusters,clust_p):
        if p_valu<=p_val:
            F_obs_plot[c]=F_obs[c]
    data=F_obs_plot.T

    #data=np.flipud(data)



        # evoked

    evok=evoked_plot(data,tmin=crop_value[0])

    fig_effect=evok.plot_image(scalings=1,units='F-value',show_names='auto',clim=dict(eeg=[0,None]),cmap='turbo')
    if png != None:
        fig_path=f'ana/results_report/images/anovas/{png}'
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        file=f'/{png}_{effect_label}_mass.png'
        filenam=fig_path+file
        txt_filename=fig_path+f'/{png}_info.txt'
        caption1=f'{file} : Time-course of {effect_label}, pval ={p_val}\n'
        with open (txt_filename,'a') as file:
            file.write(caption1)
        fig_effect.savefig(filenam,dpi=600)
    fig_topo=evok.plot_topomap(topo_times,outlines='head',scalings=1,vmin=0,cmap='turbo',units='F-value',average=plot_average)
    if png != None:
        fig_path=f'ana/results_report/images/anovas/{png}'
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        file=f'/{png}_{effect_label}_topo.png'
        filenam=fig_path+file
        txt_filename=fig_path+f'/{png}_info.txt'
        caption2=f'{file} : Topoplot of {effect_label} on {topo_times}, p val={p_val}\n'
        with open (txt_filename,'a') as file:
            file.write(caption2)
        fig_topo.savefig(filenam,dpi=600,transparent=True)


        if report != None:

            # add graphs to report to produce HTML only if label is present

            caption1=f'Time-course of {effect_label} {corr}, pval ={p_val}'
            caption2=f'Topoplot of {effect_label} on {topo_times} {corr}, p val={p_val}'

            report.add_figure(fig=fig_effect,title=f'Anova {effects_labels[0]} X {effects_labels[1]}:  {effect_label} image ',caption=caption1,image_format='svg')
            report.add_figure(fig=fig_topo,title=f'Anova {effects_labels[0]} X {effects_labels[1]}:  {effect_label} topoplot ',caption=caption2,image_format='svg')

    return report
