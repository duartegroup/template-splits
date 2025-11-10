# -----------------------------------------------------------------------------
# This file was originally distributed as part of the "LocalTemplate" project
# under the Apache License, Version 2.0 (the "License").
# 
# Copyright 2021 Shuan Chen
# 
# Modifications were made by Sara Tanovic on 10.11.2025 to adapt the code
# for the dataset and workflow used in this publication.
# 
# You may not use this file except in compliance with the License.
# A copy of the License is provided in the LICENSE file in this repository.
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------

from collections import defaultdict
import pandas as pd
import sys, os

from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')

from .template_extractor import extract_from_reaction
            
def build_template_extractor(args):
    setting = {'verbose': False, 'use_stereo': False, 'use_symbol': False, 'max_unmap': 5, 'retro': False, 'remote': True, 'least_atom_num': 2}
    for k in setting.keys():
        if k in args.keys():
            setting[k] = args[k]
    if args['retro']:
        setting['use_symbol'] = True
    print ('Template extractor setting:', setting)
    return lambda x: extract_from_reaction(x, setting)

def get_reaction_template(extractor, rxn, _id = 0):
    rxn = {'reactants': rxn.split('>>')[0], 'products': rxn.split('>>')[1], '_id': _id}
    result = extractor(rxn)
    return rxn, result

def get_full_template(template, H_change, Charge_change, Chiral_change):
    H_code = ''.join([str(H_change[k+1]) for k in range(len(H_change))])
    Charge_code = ''.join([str(Charge_change[k+1]) for k in range(len(Charge_change))])
    Chiral_code = ''.join([str(Chiral_change[k+1]) for k in range(len(Chiral_change))])
    if Chiral_code == '':
        return '_'.join([template, H_code, Charge_code])
    else:
        return '_'.join([template, H_code, Charge_code, Chiral_code])
            
def extract_templates(rxns, args, extractor):    
    TemplateEdits = {}
    TemplateCs = {}
    TemplateHs = {}
    TemplateSs = {}
    TemplateFreq = defaultdict(int)
    templates_A = defaultdict(int)
    templates_B = defaultdict(int)
    template_labels = []
    
    for i, reaction in enumerate(rxns):
        try:
            rxn, result = get_reaction_template(extractor, reaction, i)
            if 'reactants' not in result or 'reaction_smarts' not in result.keys():
                raise Exception("template problem")
            template = result['reaction_smarts']
            edits = result['edits']
            H_change = result['H_change']
            Charge_change = result['Charge_change']
            if args['use_stereo']:
                Chiral_change = result['Chiral_change']
            else:
                Chiral_change = {}
            template_H = get_full_template(template, H_change, Charge_change, Chiral_change)
            if template_H not in TemplateHs.keys():
                TemplateEdits[template_H] = {edit_type: edits[edit_type][2] for edit_type in edits}
                TemplateHs[template_H] = H_change
                TemplateCs[template_H] = Charge_change
                TemplateSs[template_H] = Chiral_change
            template_labels.append(template_H)

            TemplateFreq[template_H] += 1

            if args['retro']:
                for edit_type, bonds in edits.items():
                    bonds = bonds[0]
                    if len(bonds) > 0:
                        if edit_type in ['A', 'R']:
                            templates_A[template_H] += 1
                        else:
                            templates_B[template_H] += 1

            else:
                for edit_type, bonds in edits.items():
                    bonds = bonds[0]
                    if len(bonds) > 0:
                        if edit_type != 'A':
                            templates_A['%s_%s' % (template_H, edit_type)] += 1
                        else:
                            templates_B['%s_%s' % (template_H, edit_type)] += 1
                
        except KeyboardInterrupt:
            print('Interrupted')
            try:
                sys.exit(0)
            except SystemExit:
                os._exit(0)
        except Exception as e:
            template_labels.append("NaN")
            
        if i % 10000 == 0 and i > 0:
            print ('\r i = %s, # of template: %s, # of atom template: %s, # of bond template: %s' % (i, len(TemplateFreq), len(templates_A), len(templates_B)), end='', flush=True)
        
    print ('\n total # of template: %s' %  len(TemplateFreq))
        
    template_infos = pd.DataFrame({'Template': k, 'edit_site':TemplateEdits[k], 'change_H': TemplateHs[k], 'change_C': TemplateCs[k], 'change_S': TemplateSs[k], 'Frequency': TemplateFreq[k]} for k in TemplateHs.keys())
    
    return template_infos, template_labels
    
