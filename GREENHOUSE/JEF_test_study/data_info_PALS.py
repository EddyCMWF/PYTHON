def JULES_sources():
    PALS_dir = '/prj/GREENHOUSE/PALS_comparison/'
    # Create list of lists for JULES sources 
    #  first element of each list is the source name, 
    #  second element is a disctionary containing meta data
    sources = [ \
            ['Jvn4.5',{'data_dir':PALS_dir+'Jvn4.5/output/', \
                       'longname':'Jv4.5',                     \
                       'colour':'b',                           \
                       'gpp_name':'gpp_gb',                    \
                       'npp_name':'npp_gb',                    \
                       'resp_s_name':'resp_s_gb',              \
                       'resp_p_name':'resp_p_gb',              \
                         } ] , \
            ['Jvn4.5_notriffid',{'data_dir':PALS_dir+'Jvn4.5_notriffid/output/', \
                       'longname':'Jv4.5_notriffid',                     \
                       'colour':'c',                           \
                       'gpp_name':'gpp_gb',                    \
                       'npp_name':'npp_gb',                    \
                       'resp_s_name':'resp_s_gb',              \
                       'resp_p_name':'resp_p_gb',              \
                         } ] , \
            ['Jvn4.5_notriffid_mod',{'data_dir':PALS_dir+'Jvn4.5_notriffid/output/mod/', \
                       'longname':'Jv4.5_notriffid_mod',           \
                       'colour':'olivedrab',                   \
                       'gpp_name':'gpp_gb',                    \
                       'npp_name':'npp_gb',                    \
                       'resp_s_name':'resp_s_gb',              \
                       'resp_p_name':'resp_p_gb',              \
                         } ] , \
            ['Jvn4.5_noT_TP',{'data_dir':PALS_dir+'Jvn4.5_noT_TP/output/', \
                       'longname':'Jv4.5_noT_TP',                     \
                       'colour':'m',                           \
                       'gpp_name':'gpp_gb',                    \
                       'npp_name':'npp_gb',                    \
                       'resp_s_name':'resp_s_gb',              \
                       'resp_p_name':'resp_p_gb',              \
                         } ] , \
            ['Jvn4.5_noT_TP_lev2',{'data_dir':PALS_dir+'Jvn4.5_noT_TP_lev2/output/', \
                       'longname':'Jv4.5_noT_TP_lev2',                     \
                       'colour':'g',                           \
                       'gpp_name':'gpp_gb',                    \
                       'npp_name':'npp_gb',                    \
                       'resp_s_name':'resp_s_gb',              \
                       'resp_p_name':'resp_p_gb',              \
                         } ] , \
            ['Jvn4.3.1',{'data_dir':PALS_dir+'Jvn4.3.1/output/', \
                         'longname':'J',                     \
                         'colour':'k',                       \
                         'gpp_name':'gpp_gb',                    \
                         'npp_name':'npp_gb',                    \
                         'resp_s_name':'resp_s_gb',              \
                         'resp_p_name':'resp_p_gb',              \
                         } ] , \
            ['Jvn4.3.1-E',{'data_dir':PALS_dir+'Jvn4.3.1-E/output/', \
                             'longname':'J-E',              \
                             'colour':'saddlebrown',                     \
                             'gpp_name':'gpp_gb',                        \
                             'npp_name':'npp_nuptake_out_gb',            \
                             'resp_s_name':'co2_soil_gb',                \
                             'resp_p_name':'resp_p_gb',              \
                             } ] , \
            ['Jvn4.3.1-E-preN',{'data_dir':PALS_dir+'Jvn4.3.1-E/output/', \
                                  'longname':'J-E-preN',   \
                                  'colour':'goldenrod',                       \
                                  'gpp_name':'gpp_gb',                        \
                                  'npp_name':'npp_gb',                        \
                                  'resp_s_name':'resp_s_gb',                  \
                                  'resp_p_name':'resp_p_gb',              \
                                  } ] , \
            ['Jvn4.3.1-E-F',{'data_dir':PALS_dir+'Jvn4.3.1-E-F/output/', \
                             'longname':'J-E-F',              \
                             'colour':'b',                     \
                             'gpp_name':'gpp_gb',                        \
                             'npp_name':'npp_nuptake_out_gb',            \
                             'resp_s_name':'co2_soil_gb',                \
                             'resp_p_name':'resp_p_gb',              \
                             } ] , \
            ['Jvn4.3.1-E-F-preN',{'data_dir':PALS_dir+'Jvn4.3.1-E-F/output/', \
                                  'longname':'J-E-F-preN',   \
                                  'colour':'gold',                       \
                                  'gpp_name':'gpp_gb',                        \
                                  'npp_name':'npp_gb',                        \
                                  'resp_s_name':'resp_s_gb',                  \
                                  'resp_p_name':'resp_p_gb',              \
                                  } ] , \
            ['Jvn4.3.1-E-F_repara',{'data_dir':PALS_dir+'Jvn4.3.1-E-F_repara/output/', \
                             'longname':'J-E-F-repara',              \
                             'colour':'gold',                     \
                             'gpp_name':'gpp_gb',                        \
                             'npp_name':'npp_nuptake_out_gb',            \
                             'resp_s_name':'co2_soil_gb',                \
                             'resp_p_name':'resp_p_gb',              \
                             } ] , \
            ['Jvn4.3.1-E-F-DOC',{'data_dir':PALS_dir+'Jvn4.3.1-E-F_doc/output/', \
                                 'longname':'J-E-F-DOC',   \
                                 'colour':'cyan',                       \
                                 'gpp_name':'gpp_gb',                        \
                                 'npp_name':'npp_gb',                        \
                                 'npp_name':'npp_nuptake_out_gb',            \
                                 'resp_s_name':'co2_soil_gb',                \
                                 'resp_p_name':'resp_p_gb',              \
                                 } ] , \
            ['Jvn4.3.1-E-F-varN',{'data_dir':PALS_dir+'Jvn4.3.1-E-F_variableN/output/', \
                                 'longname':'J-E-F-varN',   \
                                 'colour':'purple',                       \
                                 'gpp_name':'gpp_gb',                        \
                                 'npp_name':'npp_gb',                        \
                                 'npp_name':'npp_nuptake_out_gb',            \
                                 'resp_s_name':'co2_soil_gb',                \
                                 'resp_p_name':'resp_p_gb',              \
                                 } ] , \
            ['Jvn4.3.1-E-F-repara-varN',{'data_dir':PALS_dir+'Jvn4.3.1-E-F_repara_variableN/output/', \
                                 'longname':'J-E-F-repara-varN',   \
                                 'colour':'orange',                       \
                                 'gpp_name':'gpp_gb',                        \
                                 'npp_name':'npp_gb',                        \
                                 'npp_name':'npp_nuptake_out_gb',            \
                                 'resp_s_name':'co2_soil_gb',                \
                                 'resp_p_name':'resp_p_gb',              \
                                 } ] , \
            ['Jvn4.3.1-E-F-lowN',{'data_dir':PALS_dir+'Jvn4.3.1-E-F_lowN/output/', \
                                 'longname':'J-E-F-lowN',   \
                                 'colour':'purple',                       \
                                 'gpp_name':'gpp_gb',                        \
                                 'npp_name':'npp_gb',                        \
                                 'npp_name':'npp_nuptake_out_gb',            \
                                 'resp_s_name':'co2_soil_gb',                \
                                 'resp_p_name':'resp_p_gb',              \
                                 } ] , \
                ]

    return sources
