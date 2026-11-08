#!/usr/bin/env python3
"""Download pinned reference artifacts and prepare a nucleotide gene index."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

RECORD = 'https://zenodo.org/api/records/21842023'


def request(url):
    if urlparse(url).scheme != 'https':
        raise ValueError('Only HTTPS URLs are supported')
    return urlopen(Request(url, headers={'User-Agent': 'RmetaG-database-helper/0.1'}), timeout=120)


def catalog():
    with request(RECORD) as response:
        record = json.load(response)
    return {item['key']: item for item in record['files'] if item['key'].endswith('.syldb')}


def fetch(url, dest, expected, source):
    if dest.exists():
        raise FileExistsError(f'Refusing to overwrite {dest}')
    manifest_path = dest.with_name(dest.name + '.manifest.json')
    if manifest_path.exists():
        raise FileExistsError(f'Refusing to overwrite {manifest_path}')
    algorithm, digest = expected.lower().split(':', 1)
    if algorithm not in ('md5', 'sha256') or len(digest) != (32 if algorithm == 'md5' else 64) or any(c not in '0123456789abcdef' for c in digest):
        raise ValueError('Expected checksum must be md5:<32 hex> or sha256:<64 hex>')
    dest.parent.mkdir(parents=True, exist_ok=True)
    checked, sha256 = hashlib.new(algorithm), hashlib.sha256()
    fd, tmp_name = tempfile.mkstemp(prefix=f'.{dest.name}.', suffix='.part', dir=dest.parent)
    temp = Path(tmp_name)
    count = 0
    try:
        with os.fdopen(fd, 'wb') as out, request(url) as response:
            while chunk := response.read(1024 * 1024):
                out.write(chunk)
                checked.update(chunk)
                sha256.update(chunk)
                count += len(chunk)
        if checked.hexdigest() != digest:
            raise ValueError(f'Checksum mismatch for {dest.name}; downloaded file discarded')
        # Refuse to replace a file created while the download was in progress.
        os.link(temp, dest)
        manifest = {'source': source, 'url': url, 'filename': dest.name, 'bytes': count,
                    'upstream_checksum': expected, 'sha256': sha256.hexdigest(),
                    'downloaded_utc': datetime.now(timezone.utc).isoformat()}
        with open(manifest_path, 'x') as handle:
            json.dump(manifest, handle, indent=2)
            handle.write('\n')
        print(f'Downloaded and verified {dest}\nManifest: {manifest_path}')
    finally:
        temp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('list-mgnify', help='List pinned MGnify-derived Sylph 0.9 indexes on Zenodo')
    get = sub.add_parser('fetch-mgnify', help='Download a biome Sylph index and verify its published checksum')
    get.add_argument('--biome', required=True)
    get.add_argument('--outdir', required=True, type=Path)
    direct = sub.add_parser('fetch-url', help='Download a reference file from HTTPS with a supplied checksum')
    direct.add_argument('--url', required=True)
    direct.add_argument('--checksum', required=True, help='sha256:<hex> or md5:<hex> from the publisher')
    direct.add_argument('--output', required=True, type=Path)
    index = sub.add_parser('index-genes', help='Build a kallisto index from nucleotide coding sequences')
    index.add_argument('--fasta', required=True, type=Path)
    index.add_argument('--output', required=True, type=Path)
    index.add_argument('--kallisto', default='kallisto', help='kallisto 0.52.0 executable')
    host = sub.add_parser('fetch-host', help='Use Hostile 2.0.2 to download a short-read host index')
    host.add_argument('--name', required=True, help='Hostile standard index name')
    host.add_argument('--outdir', required=True, type=Path)
    host.add_argument('--hostile', default='hostile')
    humann = sub.add_parser('fetch-humann', help='Use HUMAnN 3.6.1 database downloader')
    humann.add_argument('--database', choices=('chocophlan', 'uniref', 'utility_mapping'), required=True)
    humann.add_argument('--build', required=True, help='An exact build from humann_databases --available')
    humann.add_argument('--outdir', required=True, type=Path)
    humann.add_argument('--downloader', default='humann_databases')
    metaphlan = sub.add_parser('fetch-metaphlan', help='Use MetaPhlAn 4.1.1 to download an explicit index')
    metaphlan.add_argument('--index', required=True, help='Exact MetaPhlAn index name')
    metaphlan.add_argument('--outdir', required=True, type=Path)
    metaphlan.add_argument('--metaphlan', default='metaphlan')
    args = parser.parse_args()
    if args.command in ('list-mgnify', 'fetch-mgnify'):
        files = catalog()
        if args.command == 'list-mgnify':
            for key in sorted(files):
                item = files[key]
                print(f"{key[:-6]}\t{item['size']}\t{item['checksum']}")
        else:
            key = args.biome + '.syldb'
            if key not in files:
                parser.error(f'Unknown biome {args.biome}; use list-mgnify')
            entry = files[key]
            fetch(entry['links']['self'], args.outdir / key, entry['checksum'], 'doi:10.5281/zenodo.21842023; MGnify-derived Sylph 0.9.0 index')
    elif args.command == 'fetch-url':
        fetch(args.url, args.output, args.checksum, 'user-selected reference')
    elif args.command in ('fetch-host', 'fetch-humann', 'fetch-metaphlan'):
        args.outdir.mkdir(parents=True, exist_ok=True)
        if args.command == 'fetch-host':
            env = dict(os.environ, HOSTILE_CACHE_DIR=str(args.outdir.resolve()))
            cmd = [args.hostile, 'index', 'fetch', '--name', args.name, '--bowtie2']
        elif args.command == 'fetch-humann':
            env = os.environ.copy()
            cmd = [args.downloader, '--download', args.database, args.build, str(args.outdir.resolve()), '--update-config', 'no']
        else:
            env = os.environ.copy()
            cmd = [args.metaphlan, '--install', '--index', args.index, '--bowtie2db', str(args.outdir.resolve())]
        print('Running:', ' '.join(cmd), flush=True)
        subprocess.run(cmd, check=True, env=env)
        print('Download complete. Record the tool and database versions and verify compatibility before analysis.')
    else:
        if not args.fasta.is_file():
            parser.error('FASTA does not exist')
        if args.output.exists():
            parser.error('Index output already exists')
        exe = shutil.which(args.kallisto)
        if not exe:
            parser.error('kallisto executable not found')
        version = subprocess.run([exe, 'version'], check=True, text=True, capture_output=True).stdout.strip()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.rmetag-index-', dir=args.output.parent) as temporary:
            temp_index = Path(temporary) / 'genes.idx'
            subprocess.run([exe, 'index', '-i', str(temp_index), str(args.fasta.resolve())], check=True)
            os.link(temp_index, args.output)
        sha = hashlib.sha256()
        with args.fasta.open('rb') as f:
            for block in iter(lambda: f.read(1024 * 1024), b''): sha.update(block)
        manifest = {'source_fasta': str(args.fasta.resolve()), 'source_sha256': sha.hexdigest(),
                    'kallisto_version': version, 'note': 'Ensure gene IDs match the functional mapping.'}
        with open(str(args.output) + '.manifest.json', 'x') as out:
            json.dump(manifest, out, indent=2)
            out.write('\n')
        print(f'Created {args.output}')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(f'Error: {error}')
