import sys

def count_successful_assemblies(log_file):
    """
    Compte le nombre de lignes où un assemblage a réussi dans un fichier de log.
    
    Arguments:
    - log_file : str, chemin vers le fichier de log.
    
    Retourne:
    - int, nombre d'assemblages réussis.
    """
    success_count = 0

    try:
        with open(log_file, 'r') as file:
            for line in file:
                if "INFO - Assemblage" in line:
                    success_count += 1
    except FileNotFoundError:
        print(f"Erreur : Le fichier {log_file} est introuvable.")
        sys.exit(1)
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier : {e}")
        sys.exit(1)

    return success_count

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage : python count_successful_assemblies.py <log_file>")
        sys.exit(1)
    
    log_file = sys.argv[1]
    success_count = count_successful_assemblies(log_file)
    print(f"Nombre d'assemblages réussis : {success_count}")
