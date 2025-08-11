// snRNA-seq Pipeline Shiny App JavaScript

// Global variables
let pipelineRunning = false;
let logUpdateInterval = null;

// Initialize app when document is ready
$(document).ready(function() {
    console.log('snRNA-seq Pipeline App initialized');
    
    // Initialize tooltips
    $('[data-toggle="tooltip"]').tooltip();
    
    // Initialize file input styling
    initializeFileInputs();
    
    // Initialize progress tracking
    initializeProgressTracking();
    
    // Initialize notifications
    initializeNotifications();
    
    // Initialize responsive behavior
    initializeResponsiveBehavior();
});

// File input styling
function initializeFileInputs() {
    $('.file-input-wrapper').each(function() {
        const wrapper = $(this);
        const input = wrapper.find('input[type="file"]');
        const label = wrapper.find('label');
        
        input.on('change', function() {
            const fileName = this.files[0] ? this.files[0].name : 'No file selected';
            label.text(fileName);
            
            if (this.files[0]) {
                label.addClass('file-selected');
                showNotification('File selected: ' + fileName, 'success');
            } else {
                label.removeClass('file-selected');
            }
        });
    });
}

// Progress tracking
function initializeProgressTracking() {
    // Monitor pipeline progress
    $(document).on('click', '#run_pipeline', function() {
        pipelineRunning = true;
        startProgressMonitoring();
    });
    
    // Update progress bar
    function updateProgress(value) {
        $('#pipeline_progress').css('width', value + '%');
        $('#pipeline_progress').text(value + '%');
    }
}

// Start progress monitoring
function startProgressMonitoring() {
    if (logUpdateInterval) {
        clearInterval(logUpdateInterval);
    }
    
    logUpdateInterval = setInterval(function() {
        if (!pipelineRunning) {
            clearInterval(logUpdateInterval);
            return;
        }
        
        // Simulate progress updates
        const currentProgress = parseInt($('#pipeline_progress').css('width')) || 0;
        const newProgress = Math.min(currentProgress + Math.random() * 10, 90);
        updateProgress(newProgress);
    }, 2000);
}

// Notification system
function initializeNotifications() {
    // Create notification container if it doesn't exist
    if (!$('#notification-container').length) {
        $('body').append('<div id="notification-container"></div>');
    }
}

// Show notification
function showNotification(message, type = 'info', duration = 5000) {
    const container = $('#notification-container');
    const notification = $(`
        <div class="notification ${type}">
            <span class="notification-message">${message}</span>
            <button class="notification-close">&times;</button>
        </div>
    `);
    
    container.append(notification);
    
    // Auto-remove after duration
    setTimeout(function() {
        notification.fadeOut(300, function() {
            $(this).remove();
        });
    }, duration);
    
    // Manual close
    notification.find('.notification-close').on('click', function() {
        notification.fadeOut(300, function() {
            $(this).remove();
        });
    });
}

// Responsive behavior
function initializeResponsiveBehavior() {
    // Handle sidebar collapse on mobile
    $(window).on('resize', function() {
        if ($(window).width() < 768) {
            $('.sidebar').addClass('sidebar-collapse');
        } else {
            $('.sidebar').removeClass('sidebar-collapse');
        }
    });
    
    // Trigger initial check
    $(window).trigger('resize');
}

// Form validation
function validateForm() {
    let isValid = true;
    const errors = [];
    
    // Check required fields
    $('input[required], select[required]').each(function() {
        if (!$(this).val()) {
            isValid = false;
            errors.push($(this).attr('name') + ' is required');
            $(this).addClass('error');
        } else {
            $(this).removeClass('error');
        }
    });
    
    // Check file inputs
    $('input[type="file"]').each(function() {
        if ($(this).attr('required') && !$(this).get(0).files.length) {
            isValid = false;
            errors.push('Please select a file for ' + $(this).attr('name'));
            $(this).closest('.file-input-wrapper').addClass('error');
        } else {
            $(this).closest('.file-input-wrapper').removeClass('error');
        }
    });
    
    // Show errors if any
    if (!isValid) {
        showNotification('Please fix the following errors: ' + errors.join(', '), 'error', 8000);
    }
    
    return isValid;
}

// Tab navigation
function initializeTabNavigation() {
    // Store active tab in session storage
    $('a[data-toggle="tab"]').on('shown.bs.tab', function(e) {
        const tabId = $(e.target).attr('href');
        sessionStorage.setItem('activeTab', tabId);
    });
    
    // Restore active tab on page load
    const activeTab = sessionStorage.getItem('activeTab');
    if (activeTab) {
        $('a[href="' + activeTab + '"]').tab('show');
    }
}

// Data table enhancements
function initializeDataTables() {
    $('.datatable').DataTable({
        responsive: true,
        pageLength: 25,
        dom: 'Bfrtip',
        buttons: ['copy', 'csv', 'excel', 'pdf', 'print'],
        language: {
            search: "Search:",
            lengthMenu: "Show _MENU_ entries per page",
            info: "Showing _START_ to _END_ of _TOTAL_ entries",
            paginate: {
                first: "First",
                last: "Last",
                next: "Next",
                previous: "Previous"
            }
        }
    });
}

// Plot enhancements
function initializePlots() {
    // Make plots responsive
    $('.plot-container').each(function() {
        const plot = $(this);
        const resizeObserver = new ResizeObserver(function(entries) {
            entries.forEach(function(entry) {
                // Trigger plot resize if it's a plotly plot
                if (plot.find('.plotly').length) {
                    Plotly.Plots.resize(plot.find('.plotly')[0]);
                }
            });
        });
        
        resizeObserver.observe(plot[0]);
    });
}

// Keyboard shortcuts
function initializeKeyboardShortcuts() {
    $(document).on('keydown', function(e) {
        // Ctrl/Cmd + Enter to run pipeline
        if ((e.ctrlKey || e.metaKey) && e.keyCode === 13) {
            e.preventDefault();
            $('#run_pipeline').click();
        }
        
        // Ctrl/Cmd + S to save configuration
        if ((e.ctrlKey || e.metaKey) && e.keyCode === 83) {
            e.preventDefault();
            saveConfiguration();
        }
        
        // Escape to close modals
        if (e.keyCode === 27) {
            $('.modal').modal('hide');
        }
    });
}

// Configuration management
function saveConfiguration() {
    const config = {
        project: {
            name: $('#project_name').val(),
            working_dir: $('#working_dir').val()
        },
        qc: {
            min_features: parseInt($('#min_features').val()),
            min_counts: parseInt($('#min_counts').val()),
            max_features: parseInt($('#max_features').val()),
            max_counts: parseInt($('#max_counts').val()),
            max_mt_percent: parseInt($('#max_mt_percent').val())
        },
        processing: {
            normalization_method: $('#normalization_method').val(),
            n_variable_features: parseInt($('#n_variable_features').val()),
            scaling_method: $('#scaling_method').val(),
            pca_dimensions: parseInt($('#pca_dimensions').val())
        },
        clustering: {
            algorithm: $('#clustering_algorithm').val(),
            resolution: parseFloat($('#clustering_resolution').val()),
            min_cluster_size: parseInt($('#min_cluster_size').val())
        }
    };
    
    // Save to localStorage
    localStorage.setItem('pipeline_config', JSON.stringify(config));
    showNotification('Configuration saved locally', 'success');
}

function loadConfiguration() {
    const savedConfig = localStorage.getItem('pipeline_config');
    if (savedConfig) {
        const config = JSON.parse(savedConfig);
        
        // Populate form fields
        $('#project_name').val(config.project.name);
        $('#working_dir').val(config.project.working_dir);
        $('#min_features').val(config.qc.min_features);
        $('#min_counts').val(config.qc.min_counts);
        $('#max_features').val(config.qc.max_features);
        $('#max_counts').val(config.qc.max_counts);
        $('#max_mt_percent').val(config.qc.max_mt_percent);
        $('#normalization_method').val(config.processing.normalization_method);
        $('#n_variable_features').val(config.processing.n_variable_features);
        $('#scaling_method').val(config.processing.scaling_method);
        $('#pca_dimensions').val(config.processing.pca_dimensions);
        $('#clustering_algorithm').val(config.clustering.algorithm);
        $('#clustering_resolution').val(config.clustering.resolution);
        $('#min_cluster_size').val(config.clustering.min_cluster_size);
        
        showNotification('Configuration loaded from local storage', 'info');
    }
}

// Export functionality
function exportResults() {
    // Create zip file with results
    const zip = new JSZip();
    
    // Add configuration
    const config = {
        timestamp: new Date().toISOString(),
        parameters: {
            project_name: $('#project_name').val(),
            min_features: $('#min_features').val(),
            min_counts: $('#min_counts').val(),
            // ... other parameters
        }
    };
    
    zip.file('configuration.json', JSON.stringify(config, null, 2));
    
    // Generate and download zip
    zip.generateAsync({type: 'blob'}).then(function(content) {
        const link = document.createElement('a');
        link.href = URL.createObjectURL(content);
        link.download = 'pipeline_results.zip';
        link.click();
    });
}

// Error handling
function handleError(error, context = '') {
    console.error('Error in ' + context + ':', error);
    showNotification('An error occurred: ' + error.message, 'error', 10000);
}

// Performance monitoring
function startPerformanceMonitoring() {
    const startTime = performance.now();
    
    return {
        end: function() {
            const endTime = performance.now();
            const duration = endTime - startTime;
            console.log('Operation completed in ' + duration.toFixed(2) + 'ms');
            return duration;
        }
    };
}

// Initialize all components
$(document).ready(function() {
    initializeTabNavigation();
    initializeDataTables();
    initializePlots();
    initializeKeyboardShortcuts();
    
    // Load saved configuration
    loadConfiguration();
    
    // Auto-save configuration on form changes
    $('input, select').on('change', function() {
        setTimeout(saveConfiguration, 1000);
    });
    
    console.log('All components initialized successfully');
});
