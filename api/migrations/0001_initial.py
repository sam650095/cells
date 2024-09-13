# Generated by Django 4.2.13 on 2024-09-13 06:26

from django.db import migrations, models


class Migration(migrations.Migration):
    initial = True

    dependencies = []

    operations = [
        migrations.CreateModel(
            name="OperationStep",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("session_id", models.PositiveIntegerField()),
                ("step", models.CharField(max_length=100)),
                ("operation_type", models.CharField(max_length=100)),
                ("input_values", models.JSONField()),
                ("output_values", models.JSONField(blank=True, null=True)),
            ],
        ),
    ]
