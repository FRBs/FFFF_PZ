# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

YSE_PZ (Young Supernova Experiment - Photometric Redshift) is a Django-based transient survey management platform that has been extended for Fast Radio Burst (FRB) follow-up. This repository appears to be a fork focused on FRB research (FFFF-PZ), building on the original YSE_PZ infrastructure. The platform ingests transient discovery alerts, identifies host galaxies, downloads archival data, and provides tools for follow-up decision making.

## Development Setup

### Docker (Recommended)

The strongly recommended way to run this project is using Docker:

1. **Create settings file**: Copy `YSE_PZ/public_settings.ini` to `YSE_PZ/settings.ini`
   ```bash
   cp YSE_PZ/public_settings.ini YSE_PZ/settings.ini
   ```

2. **Create environment file**: Copy `docker/public.env` to `docker/.env`
   ```bash
   cp docker/public.env docker/.env
   ```
   Edit `.env` to set required variables: `VOL`, `VOL_DB`, `DB_PWD`, `STATIC_VOL`, `LOCAL_DB_PORT`, `LOCAL_HTTP_PORT`, `DJANGO_SUPERUSER_USERNAME`, `DJANGO_SUPERUSER_PASSWORD`, `DJANGO_SUPERUSER_EMAIL`

3. **Start containers**:
   ```bash
   cd docker
   docker compose up
   ```

4. **Collect static files** (first run):
   ```bash
   docker exec -it ysepz_web_container bash
   python3 manage.py collectstatic
   exit
   ```

5. **Create superuser** (if not auto-created):
   ```bash
   docker exec -it ysepz_web_container bash -c 'python3 manage.py createsuperuser --noinput'
   ```

6. **Access the application**: Navigate to `http://0.0.0.0/login/`

7. **Stop containers**:
   ```bash
   cd docker
   docker compose down
   ```

### Native Installation (Not Recommended)

If needed, native installation requires:
- MySQL 8.0+
- Anaconda Python environment
- See `docs/install.rst` for detailed instructions

## Common Commands

### Django Management

- **Run development server** (native): `python manage.py runserver`
- **Apply migrations**: `python manage.py migrate`
- **Create migrations**: `python manage.py makemigrations`
- **Create superuser**: `python manage.py createsuperuser`
- **Collect static files**: `python manage.py collectstatic`
- **Django shell**: `python manage.py shell`
- **Database shell**: `python manage.py dbshell`

### Running Cron Jobs (Data Ingest)

The platform uses django-cron for data ingestion. To manually run cron jobs:

```bash
# Ingest recent TNS transients
python manage.py runcrons YSE_App.data_ingest.TNS_uploads.TNS_recent --force

# Update TNS data for followed transients
python manage.py runcrons YSE_App.data_ingest.TNS_uploads.TNS_updates --force

# Update ignored transients
python manage.py runcrons YSE_App.data_ingest.TNS_uploads.TNS_Ignore_updates --force
```

All available cron classes are listed in `YSE_PZ/settings.py` under `CRON_CLASSES`.

### Docker-Specific Commands

- **Enter web container**: `docker exec -it ysepz_web_container bash`
- **Enter database container**: `docker exec -it ysepz_db_container bash`
- **View logs**: `docker compose logs -f`
- **Restart services**: `docker compose restart`

## Architecture Overview

### Core Django Structure

This is a Django 1.11+ application with the following structure:

- **`YSE_PZ/`**: Main Django project directory containing settings, URLs, and WSGI configuration
  - `settings.py`: Loads configuration from `settings.ini` file
  - `urls.py`: Root URL configuration (delegates to YSE_App)

- **`YSE_App/`**: Primary Django application containing all models, views, and business logic
  - Entry point: `manage.py` (standard Django management script)

### Key Application Components

#### Models (`YSE_App/models/`)

Models are organized by domain into separate files:

- **Core transient models**: `transient_models.py`, `host_models.py`, `tag_models.py`
- **FRB-specific models**: `frbtransient_models.py`, `frbgalaxy_models.py`, `frbphot_models.py`, `frbfollowup_models.py`, `frbsample_models.py`
- **Observation models**: `followup_models.py`, `observation_task_models.py`, `telescope_models.py`, `instrument_models.py`
- **Data models**: `phot_models.py`, `spectra_models.py`, `photometric_band_models.py`
- **System models**: `enum_models.py`, `log_models.py`, `profile_models.py`
- **Survey models**: `survey_models.py`, `gw_models.py`

All models are imported through `models/__init__.py`.

#### Views and APIs

- **`views.py`**: Main web interface views (dashboard, transient detail pages)
- **`yse_views.py`**: YSE-specific views for supernovae
- **`frb_*.py`** files: FRB-specific functionality
  - `frb_init.py`: Database initialization for FRB objects
  - `frb_status.py`: FRB status management
  - `frb_tags.py`: FRB tagging system
  - `frb_targeting.py`: Target selection logic
  - `frb_observing.py`: Observing plan generation
  - `frb_utils.py`: FRB utility functions
  - `frb_tables.py`: Table rendering for FRB data
- **`api_views.py`**: REST API viewsets using Django REST Framework
- **`form_views.py`**: Form-based views for data entry
- **`data_utils.py`**: Data manipulation and utilities
- **`urls.py`**: URL routing (includes FRB-specific routes under `/frb_*` paths)

#### Data Ingestion (`YSE_App/data_ingest/`)

Automated data ingestion through cron jobs:

- **`TNS_uploads.py`**: Transient Name Server integration (creates/updates transients)
- **`Query_ZTF.py`**: ZTF/MARS photometry ingestion
- **`QUB_data.py`**: Queen's University Belfast data ingestion
- **`Photo_Z.py`**, **`PS1_PhotoZ.py`**, **`SDSS_Photo_Z.py`**: Photometric redshift calculations
- **`host_associate.py`**: Host galaxy association
- **`YSE_Forced_Phot.py`**: Forced photometry processing
- **`DECam_upload.py`**: DECam data uploads
- **`Gaia_LC.py`**: Gaia light curve ingestion

#### Serializers (`YSE_App/serializers/`)

Django REST Framework serializers for API endpoints, organized by model type (transient, followup, instrument, etc.).

#### Templates (`YSE_App/templates/`)

Django templates for web interface rendering.

### FRB-Specific Architecture

The codebase has been extended with FRB (Fast Radio Burst) functionality:

1. **FRB Models**: Parallel model hierarchy to standard transients
   - `FRBTransient`: Main FRB object with position, DM (Dispersion Measure), etc.
   - `FRBGalaxy`: Host galaxy candidates
   - `FRBFollowupRequest`: Follow-up observation requests
   - `FRBSampleCriteria`: Sample selection criteria

2. **FRB Views**: Dedicated dashboard and detail views (`/frb_dashboard/`, `/frb_transient_detail/<name>/`)

3. **FRB Utilities**: Helper modules (`frb_utils.py`, `frb_init.py`, `frb_status.py`, `frb_tags.py`) for FRB-specific operations

4. **CHIME Integration**: Special handling for CHIME/FRB survey data (`YSE_App/chime/`)

### Database Configuration

- Uses MySQL 8.0 as the database backend
- Multiple database connections defined:
  - `default`: Main application database
  - `explorer`: Read-only database for SQL Explorer queries
- Configuration loaded from `settings.ini` file
- Database migrations tracked in `YSE_App/migrations/`

### Authentication and Permissions

- Uses Django's built-in authentication system
- REST API uses token authentication and session authentication
- Group-based permissions for multi-collaboration support
- Audit logging enabled via `django-auditlog`

### Static Files and Frontend

- Static files served through nginx in Docker setup
- Uses Bootstrap, jQuery, Chart.js, and various other frontend libraries
- Frontend assets located in `YSE_App/static/YSE_App/`

## Important Conventions

### Database Operations

- Use `frb_utils.add_or_grab_obj()` for FRB object creation to avoid duplicates
- Models inherit from `BaseModel` which provides common fields (created_at, modified_at, created_by, modified_by)
- Foreign keys use `models.SET()` with sentinel functions to handle deletions

### FRB Naming

- FRB names follow the format: `FRB<YYYYMMDD><letter>` (e.g., `FRB20220912A`)
- Names must be unique in the database

### Status Management

- FRB transients have status workflow defined in `frb_status.py`
- Status changes are tracked through `TransientStatus` model

### API Design

- API endpoints follow REST conventions
- ViewSets extend `custom_viewsets.ListCreateRetrieveUpdateViewSet` from `YSE_App/common/`
- Filtering provided through `django-filters`

## Configuration Notes

### Settings File (`YSE_PZ/settings.ini`)

Critical configuration includes:
- Database credentials and connection info
- SMTP settings for email notifications
- External service credentials (LCOGT, TNS, ZTF, etc.)
- Site-specific settings (DEBUG, STATIC_URL, etc.)

### Environment Variables (Docker)

The Docker setup uses `.env` file for:
- Volume mappings
- Port assignments
- Database passwords
- Superuser credentials

## Testing

- Test directory: `YSE_App/tests/`
- Currently minimal test coverage
- Run tests with: `python manage.py test YSE_App`

## Common Pitfalls

1. **Settings file location**: Must be at `YSE_PZ/settings.ini`, not in the repository root
2. **Docker volumes**: Ensure `VOL` paths in `.env` are absolute paths
3. **Static files**: Must run `collectstatic` after code changes affecting static assets
4. **Database migrations**: After model changes, always create and apply migrations
5. **Host validation**: Custom regex in `manage.py` relaxes host validation for development

## Documentation

- Full documentation: https://yse-pz.readthedocs.io/
- Documentation source: `docs/` directory (Sphinx/RST format)
- Build docs: `cd docs && make html`
