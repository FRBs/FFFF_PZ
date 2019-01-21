# Generated by Django 2.0.4 on 2019-01-19 00:59

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('explorer', '0008_auto_20180911_2051'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('YSE_App', '0023_photometricband_disp_symbol'),
    ]

    operations = [
        migrations.CreateModel(
            name='UserQuery',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_date', models.DateTimeField(auto_now_add=True)),
                ('modified_date', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(max_length=64)),
                ('created_by', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, related_name='userquery_created_by', to=settings.AUTH_USER_MODEL)),
                ('modified_by', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, related_name='userquery_modified_by', to=settings.AUTH_USER_MODEL)),
                ('query', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='explorer.Query')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'abstract': False,
            },
        ),
    ]