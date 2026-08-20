#!/usr/bin/env python3
import os
import glob
import json
import warnings
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from Bio.PDB import MMCIFParser, NeighborSearch

RESIDUOS_PROTEINA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
RESIDUOS_ARN = {'A', 'U', 'G', 'C', 'RA', 'RU', 'RG', 'RC', 'DA', 'DT', 'DG', 'DC'}

def calcular_interacciones_cif(ruta_cif, dist_min, dist_max):
    parser = MMCIFParser(QUIET=True)
    estructura = parser.get_structure('modelo', ruta_cif)
    
    atomos_proteina, atomos_arn = [], []
    for modelo in estructura:
        for cadena in modelo:
            for residuo in cadena:
                nombre_res = residuo.get_resname().strip().upper()
                if nombre_res in RESIDUOS_PROTEINA:
                    atomos_proteina.extend(residuo.get_atoms())
                elif nombre_res in RESIDUOS_ARN:
                    atomos_arn.extend(residuo.get_atoms())
                    
    if not atomos_proteina or not atomos_arn:
        return 0, 0, {}

    motor_busqueda = NeighborSearch(atomos_proteina)
    distancias_pares = {}
    choques_estericos = set()
    
    for atomo_arn in atomos_arn:
        atomos_cercanos = motor_busqueda.search(atomo_arn.get_coord(), dist_max)
        for atomo_prot in atomos_cercanos:
            distancia_real = atomo_arn - atomo_prot
            res_prot, res_arn = atomo_prot.get_parent(), atomo_arn.get_parent()
            
            id_prot = f"{res_prot.get_resname().strip()} {res_prot.get_id()[1]}"
            id_arn = f"{res_arn.get_resname().strip()} {res_arn.get_id()[1]}"
            pareja = (id_prot, id_arn)
            
            if distancia_real < dist_min:
                choques_estericos.add(pareja)
            else:
                if pareja not in distancias_pares or distancia_real < distancias_pares[pareja]:
                    distancias_pares[pareja] = distancia_real
                
    interacciones_validas = {pareja: dist for pareja, dist in distancias_pares.items() if pareja not in choques_estericos}
    return len(interacciones_validas), len(choques_estericos), interacciones_validas

def ejecutar_pipeline(carpeta_raiz, umbral_iptm, dist_min, dist_max, nombre_proteina):
    print(f"Análisis en curso")
    print(f"Proteína objetivo: {nombre_proteina}")
    print(f"Analizando directorio: {carpeta_raiz}")
    print(f"Umbral ipTM configurado: {umbral_iptm}")
    print(f"Distancias de Interacción: {dist_min} Å (Mín) - {dist_max} Å (Máx)\n")
    
    # --- 1: Índice General ---
    print("1. Extrayendo datos de archivos JSON...")
    datos_indice = []
    archivos_json = glob.glob(os.path.join(carpeta_raiz, '**', '*.json'), recursive=True)
    
    for ruta_json in archivos_json:
        try:
            with open(ruta_json, 'r') as f:
                dj = json.load(f)
                
            ptm, iptm = dj.get('ptm'), dj.get('iptm')
            if ptm is None or iptm is None: continue
                
            dir_json = os.path.dirname(ruta_json)
            cifs = glob.glob(os.path.join(dir_json, '*.cif'))
            if not cifs: continue
                
            datos_indice.append({
                'Carpeta': os.path.basename(dir_json),
                'Nombre_CIF': os.path.basename(cifs[0]),
                'Ruta_CIF': cifs[0],
                'pTM': float(ptm),
                'ipTM': float(iptm)
            })
        except Exception as e:
            pass

    df_indice = pd.DataFrame(datos_indice)
    if df_indice.empty:
        print("No se encontraron modelos válidos.")
        return
        
    df_indice.to_csv(os.path.join(carpeta_raiz, f"01_indice_general_modelos_{nombre_proteina}.csv"), index=False)

    # --- Histogramas Iniciales ---
    print("2. Generando histogramas de confianza")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    ax1.hist(df_indice['pTM'], bins=40, color='#4C72B0', edgecolor='black', alpha=0.8)
    ax1.set_title('Distribución Global de pTM', weight='bold')
    ax1.set_xticks(np.arange(0, 1.05, 0.05))
    ax1.tick_params(axis='x', rotation=45)
    
    ax2.hist(df_indice['ipTM'], bins=40, color='#55A868', edgecolor='black', alpha=0.8)
    ax2.axvline(umbral_iptm, color='red', linestyle='dashed', linewidth=2, label=f'Corte: {umbral_iptm}')
    ax2.set_title('Distribución Global de ipTM', weight='bold')
    ax2.set_xticks(np.arange(0, 1.05, 0.05))
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(carpeta_raiz, f"02_histograma_inicial_ptm_iptm_{nombre_proteina}.png"), dpi=300)
    plt.close()

    # --- Filtrado ipTM ---
    print(f"3. Corte de ipTM >= {umbral_iptm}...")
    df_filtrado = df_indice[df_indice['ipTM'] >= umbral_iptm].copy()
    print(f"Modelos filtrados: {len(df_filtrado)}")
    
    if df_filtrado.empty: return
    
    nombre_csv_filtrado = f"03_indice_modelos_filtrados_{nombre_proteina}_ipTM{umbral_iptm}.csv"
    nombre_png_filtrado = f"04_histograma_filtrado_{nombre_proteina}_ipTM{umbral_iptm}.png"
        
    df_filtrado.to_csv(os.path.join(carpeta_raiz, nombre_csv_filtrado), index=False)

    plt.figure(figsize=(8, 5))
    plt.hist(df_filtrado['ipTM'], bins=30, color='#DD8452', edgecolor='black')
    plt.title(f'Distribución ipTM ({nombre_proteina} | Modelos >= {umbral_iptm})', weight='bold')
    plt.xticks(np.arange(umbral_iptm, 1.05, 0.02), rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(carpeta_raiz, nombre_png_filtrado), dpi=300)
    plt.close()

    # --- Análisis Estructural ---
    print("4. Analizando distancias")
    res_int, det_enl = [], []
    
    for idx, row in df_filtrado.iterrows():
        n_val, n_choq, dic_int = calcular_interacciones_cif(row['Ruta_CIF'], dist_min, dist_max)
        
        res_int.append({**row.to_dict(), 'Interacciones_Validas': n_val, 'Choques_Descartados': n_choq})
        for (r_prot, r_arn), dist in dic_int.items():
            det_enl.append({'Carpeta': row['Carpeta'], 'Nombre_CIF': row['Nombre_CIF'], 
                            'Residuo_Proteina': r_prot, 'Residuo_ARN': r_arn, 'Distancia_A': round(dist, 3)})

    df_resumen_int = pd.DataFrame(res_int)
    df_resumen_int.to_csv(os.path.join(carpeta_raiz, f"05_resumen_interacciones_{nombre_proteina}.csv"), index=False)
    
    df_detalles_int = pd.DataFrame(det_enl)
    if not df_detalles_int.empty:
        df_detalles_int.to_csv(os.path.join(carpeta_raiz, f"06_detalles_distancias_{nombre_proteina}.csv"), index=False)

    # --- Histograma Final Interacciones ---
    print("5. Generando histograma de distancias")
    plt.figure(figsize=(10, 6))
    max_int = df_resumen_int['Interacciones_Validas'].max()
    
    plt.hist(df_resumen_int['Interacciones_Validas'], bins=range(0, (max_int or 10) + 2), color='#8172B3', edgecolor='black', align='left')
    plt.title(f'Interacciones Válidas ({nombre_proteina} | {dist_min}Å - {dist_max}Å)', weight='bold')
    plt.savefig(os.path.join(carpeta_raiz, f"07_histograma_interacciones_{nombre_proteina}.png"), dpi=300)
    plt.close()

    # --- Resumen de Residuos Críticos ---
    print("6. Calculando distancias mínimas por residuo específico")
    if not df_detalles_int.empty:
        df_prot = df_detalles_int.groupby('Residuo_Proteina').agg(
            Apariciones=('Residuo_Proteina', 'count'),
            Distancia_Minima=('Distancia_A', 'min')
        ).reset_index().sort_values('Distancia_Minima')
        df_prot.to_csv(os.path.join(carpeta_raiz, f"08_residuos_proteina_min_distancia_{nombre_proteina}.csv"), index=False)

        df_arn = df_detalles_int.groupby('Residuo_ARN').agg(
            Apariciones=('Residuo_ARN', 'count'),
            Distancia_Minima=('Distancia_A', 'min')
        ).reset_index().sort_values('Distancia_Minima')
        df_arn.to_csv(os.path.join(carpeta_raiz, f"09_residuos_arn_min_distancia_{nombre_proteina}.csv"), index=False)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        ax1.hist(df_prot['Distancia_Minima'], bins=15, color='#D55E00', edgecolor='black', alpha=0.8)
        ax1.set_title(f'Distancias Mínimas ({nombre_proteina})', weight='bold')
        ax1.set_xlabel('Distancia Mínima de Interacción (Å)')
        ax1.set_ylabel('Frecuencia')
        ax1.grid(axis='y', linestyle='--', alpha=0.7)

        ax2.hist(df_arn['Distancia_Minima'], bins=15, color='#0072B2', edgecolor='black', alpha=0.8)
        ax2.set_title(f'Distancias Mínimas (ARN asociado a {nombre_proteina})', weight='bold')
        ax2.set_xlabel('Distancia Mínima de Interacción (Å)')
        ax2.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        plt.savefig(os.path.join(carpeta_raiz, f"10_histogramas_distancias_minimas_residuos_{nombre_proteina}.png"), dpi=300)
        plt.close()

    print("\n=== PIPELINE FINALIZADO CORRECTAMENTE ===")
    print(f"Los resultados se han guardado dentro de: {os.path.abspath(carpeta_raiz)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline de Análisis Blitz para complejos Proteína-ARN")
    
    parser.add_argument("-i", "--input", required=True, help="Ruta a la carpeta de simulaciones.")
    parser.add_argument("-u", "--umbral", type=float, default=0.75, help="Corte para ipTM (default: 0.75)")
    parser.add_argument("-p", "--proteina", type=str, default="Proteina", help="Nombre de la proteína para etiquetar archivos")
    parser.add_argument("--dmin", type=float, default=2.0, help="Distancia mínima en Å para descartar choques (default: 2.0)")
    parser.add_argument("--dmax", type=float, default=3.5, help="Distancia máxima en Å para interacción (default: 3.5)")
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input):
        print(f"Error: El directorio '{args.input}' no existe.")
    else:
        ejecutar_pipeline(args.input, args.umbral, args.dmin, args.dmax, args.proteina)
